# gpharm_core/named_allele_matcher.py

import json
import logging
import itertools
from pathlib import Path
from dataclasses import dataclass, field

logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# Типы данных
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class VariantCall:
    """Один вариант из VCF."""
    rsid: str
    chrom: str
    pos: int
    ref: str
    alt: list[str]       # список альтернативных аллелей
    genotype: tuple[str, str]  # например ('C', 'T') — всегда 2 аллеля
    phased: bool         # True если GT разделён '|', False если '/'
    quality: str         # FILTER поле


@dataclass
class NamedAlleleDefinition:
    """Один star-аллель из CPIC definition table."""
    name: str                          # например '*4'
    defining_variants: dict[str, str]  # {rsid: expected_base}  — только непустые ячейки
    score: int = 0                     # число defining_variants (заполняется при загрузке)


@dataclass
class DiplotypeCall:
    """Финальный результат для одного гена — контракт с Dev 2."""
    gene: str
    allele_a: str                       # например '*2'
    allele_b: str                       # например '*4'
    phased: bool
    confidence: float                   # 0.0 – 1.0
    diplotype_score: int                # сумма score(a) + score(b)
    missing_positions: list[str] = field(default_factory=list)  # rsID без данных
    low_quality_positions: list[str] = field(default_factory=list)
    all_candidates: list[tuple[str, str, int]] = field(default_factory=list)  # (a, b, score)


# ─────────────────────────────────────────────────────────────────────────────
# Шаг 1: парсинг VCF → {rsID: VariantCall}
# ─────────────────────────────────────────────────────────────────────────────

def extract_pgx_variants(vcf_path: str | Path) -> dict[str, VariantCall]:
    """
    Парсит VCF-файл и возвращает словарь {rsID: VariantCall}.

    Что делает:
    - Читает VCF построчно, пропускает заголовки (#)
    - Для каждого варианта восстанавливает реальные нуклеотиды из GT-индексов
    - Обрабатывает фазированные ('|') и нефазированные ('/') генотипы
    - Пропускает варианты без rsID (ID == '.')

    Аргументы:
        vcf_path: путь к VCF (plain, после preprocess_vcf)

    Возвращает:
        dict {rsid: VariantCall} — только варианты с непустым rsID
    """
    vcf_path = Path(vcf_path)
    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF не найден: {vcf_path}")

    variants: dict[str, VariantCall] = {}
    format_col_idx: int | None = None  # индекс GT в FORMAT

    with open(vcf_path, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")

            # Пропускаем заголовки
            if line.startswith("##"):
                continue

            # Строка с именами колонок — определяем индекс GT в FORMAT
            if line.startswith("#CHROM"):
                cols = line.lstrip("#").split("\t")
                # Ожидаем стандартный VCF: CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE
                continue

            cols = line.split("\t")
            if len(cols) < 10:
                continue  # некорректная строка

            chrom   = cols[0]
            pos     = int(cols[1])
            rsid    = cols[2]
            ref     = cols[3]
            alt_raw = cols[4]      # может быть 'A,T' для multi-allelic (уже split после norm)
            quality = cols[6]      # FILTER
            fmt     = cols[8]      # FORMAT: 'GT:DP:GQ...'
            sample  = cols[9]      # SAMPLE значения

            # Пропускаем варианты без rsID
            if rsid == "." or not rsid.startswith("rs"):
                continue

            # Список ALT-аллелей
            alt_alleles = alt_raw.split(",")

            # Индекс GT в FORMAT (кешируем)
            fmt_fields = fmt.split(":")
            try:
                gt_idx = fmt_fields.index("GT")
            except ValueError:
                logger.warning(f"Нет поля GT в FORMAT для {rsid}, пропускаем")
                continue

            # Значение GT из SAMPLE
            sample_fields = sample.split(":")
            if gt_idx >= len(sample_fields):
                continue
            gt_raw = sample_fields[gt_idx]  # например '0/1', '1|0', './.'

            # Определяем фазирование
            phased = "|" in gt_raw
            sep = "|" if phased else "/"
            gt_parts = gt_raw.split(sep)

            # Пропускаем missing (./., .|.)
            if "." in gt_parts:
                logger.debug(f"Missing genotype для {rsid}: {gt_raw}")
                continue

            # Конвертируем индексы (0=REF, 1=ALT[0], 2=ALT[1], ...) в нуклеотиды
            allele_map = [ref] + alt_alleles  # индекс → нуклеотид
            try:
                allele_a = allele_map[int(gt_parts[0])]
                allele_b = allele_map[int(gt_parts[1])]
            except (IndexError, ValueError) as e:
                logger.warning(f"Не удалось разобрать GT '{gt_raw}' для {rsid}: {e}")
                continue

            variants[rsid] = VariantCall(
                rsid=rsid,
                chrom=chrom,
                pos=pos,
                ref=ref,
                alt=alt_alleles,
                genotype=(allele_a, allele_b),
                phased=phased,
                quality=quality,
            )

    logger.info(f"Извлечено {len(variants)} вариантов с rsID из {vcf_path.name}")
    return variants


# ─────────────────────────────────────────────────────────────────────────────
# Вспомогательные функции для Named Allele Matcher
# ─────────────────────────────────────────────────────────────────────────────

def _load_allele_definitions(json_path: Path) -> tuple[str, list[NamedAlleleDefinition]]:
    """
    Читает CPIC allele definition JSON для одного гена.

    Формат CPIC JSON (PharmCAT):
        {
          "gene": "CYP2C19",
          "variants": [{"rsid": "rs4244285", ...}, ...],
          "namedAlleles": [
            {"name": "*1", "alleles": ["G", "C", null, ...]},
            {"name": "*2", "alleles": ["A", null, null, ...]}
          ]
        }

    null в alleles означает «любой» (позиция не определяет этот аллель),
    реальное значение — конкретный нуклеотид, который ДОЛЖЕН быть у пациента.
    """
    with open(json_path, "r") as f:
        data = json.load(f)

    gene = data["gene"]

    # Список rsID в том же порядке что alleles[]
    variant_rsids: list[str] = [
        v.get("rsid", v.get("chromosomeHgvsName", f"pos_{i}"))
        for i, v in enumerate(data["variants"])
    ]

    named_alleles: list[NamedAlleleDefinition] = []
    for na in data["namedAlleles"]:
        defining: dict[str, str] = {}
        for rsid, base in zip(variant_rsids, na["alleles"]):
            if base is not None:  # null → позиция не задействована для этого аллеля
                defining[rsid] = base
        allele = NamedAlleleDefinition(
            name=na["name"],
            defining_variants=defining,
            score=len(defining),
        )
        named_alleles.append(allele)

    logger.debug(f"Загружено {len(named_alleles)} аллелей для {gene}")
    return gene, named_alleles


def _match_single_haplotype(
    haplotype: dict[str, str],        # {rsid: нуклеотид} — одна хромосома
    allele_def: NamedAlleleDefinition,
    available_rsids: set[str],        # позиции которые вообще есть в VCF
) -> bool:
    """
    Проверяет: совпадает ли хаплотип пациента с определением аллеля?

    Правило CPIC:
    - Если позиция определяет аллель (не null) И есть в VCF → должна совпасть
    - Если позиция определяет аллель НО отсутствует в VCF → пропускаем (ref assumed)
    - Позиции с null не проверяются вообще
    """
    for rsid, expected_base in allele_def.defining_variants.items():
        if rsid not in available_rsids:
            continue  # missing позиция → считаем ref, не штрафуем
        observed = haplotype.get(rsid)
        if observed is None:
            continue
        if observed != expected_base:
            return False  # явное несовпадение
    return True


def _compute_allele_score(
    allele_def: NamedAlleleDefinition,
    available_rsids: set[str],
) -> int:
    """
    Скор = число defining positions, которые реально присутствуют в VCF.
    Именно так работает PharmCAT scoring: missing позиции снижают скор.
    """
    return sum(1 for rsid in allele_def.defining_variants if rsid in available_rsids)


# ─────────────────────────────────────────────────────────────────────────────
# Шаг 2: Named Allele Matcher — главная логика
# ─────────────────────────────────────────────────────────────────────────────

def call_star_alleles(
    variants: dict[str, VariantCall],
    allele_def_dir: str | Path,
    gene: str | None = None,
) -> dict[str, DiplotypeCall]:
    """
    Named Allele Matcher: для каждого гена определяет наиболее вероятный диплотип.

    Алгоритм (точное воспроизведение логики PharmCAT):

    [ФАЗИРОВАННЫЕ данные]
        1. Для каждой хромосомы (A и B) берём её нуклеотид по каждому rsID
        2. Перебираем все named alleles, проверяем _match_single_haplotype()
        3. Из совпавших выбираем с максимальным скором
        4. Пара (best_A, best_B) = диплотип

    [НЕФАЗИРОВАННЫЕ данные]
        1. Для каждого rsID у нас есть пара нуклеотидов без информации о хромосоме
        2. Генерируем все перестановки распределения аллелей по хромосомам
           (2^N комбинаций для N гетерозиготных позиций — ограничено до 2^12)
        3. Для каждой перестановки ищем совместимые диплотипы
        4. Проверяем consistency: аллели двух хромосом должны «объяснить» все варианты
        5. Выбираем диплотип с максимальным суммарным скором

    Аргументы:
        variants:       результат extract_pgx_variants()
        allele_def_dir: папка с CPIC JSON файлами (CYP2C19.json, CYP2D6.json, ...)
        gene:           если задан — обрабатывает только этот ген; иначе все JSON в папке

    Возвращает:
        dict {gene_name: DiplotypeCall}
    """
    allele_def_dir = Path(allele_def_dir)
    available_rsids = set(variants.keys())

    # Определяем какие файлы обрабатывать
    if gene:
        json_files = [allele_def_dir / f"{gene}.json"]
    else:
        json_files = sorted(allele_def_dir.glob("*.json"))

    results: dict[str, DiplotypeCall] = {}

    for json_path in json_files:
        if not json_path.exists():
            logger.warning(f"Файл определений не найден: {json_path}")
            continue

        gene_name, named_alleles = _load_allele_definitions(json_path)

        # rsID, которые нужны для этого гена
        gene_rsids: set[str] = set()
        for na in named_alleles:
            gene_rsids.update(na.defining_variants.keys())

        # Позиции, нужные гену, но отсутствующие в VCF
        missing_positions = sorted(gene_rsids - available_rsids)
        gene_available = gene_rsids & available_rsids
        low_qual = [
            r for r in gene_available
            if variants[r].quality not in ("PASS", ".", "")
        ]

        if not gene_available:
            logger.warning(f"{gene_name}: нет ни одной позиции в VCF → no-call")
            continue

        logger.info(
            f"{gene_name}: {len(gene_available)}/{len(gene_rsids)} позиций покрыто "
            f"({len(missing_positions)} missing)"
        )

        # Определяем фазирование: если ВСЕ нужные варианты фазированы — используем phased
        is_phased = all(
            variants[r].phased for r in gene_available if r in variants
        )

        diplotype: DiplotypeCall

        if is_phased:
            diplotype = _call_phased(
                gene_name, named_alleles, variants, gene_available, missing_positions, low_qual
            )
        else:
            diplotype = _call_unphased(
                gene_name, named_alleles, variants, gene_available, missing_positions, low_qual
            )

        results[gene_name] = diplotype
        logger.info(
            f"{gene_name}: {diplotype.allele_a}/{diplotype.allele_b} "
            f"(score={diplotype.diplotype_score}, conf={diplotype.confidence:.2f})"
        )

    return results


# ─────────────────────────────────────────────────────────────────────────────
# Внутренние функции вызова
# ─────────────────────────────────────────────────────────────────────────────

def _call_phased(
    gene: str,
    named_alleles: list[NamedAlleleDefinition],
    variants: dict[str, VariantCall],
    gene_available: set[str],
    missing_positions: list[str],
    low_qual: list[str],
) -> DiplotypeCall:
    """
    Вызов диплотипа для ФАЗИРОВАННЫХ данных.
    Каждая хромосома матчится независимо.
    """
    # Строим словари хаплотипов: {rsid: нуклеотид} для каждой хромосомы
    hap_a: dict[str, str] = {}
    hap_b: dict[str, str] = {}

    for rsid in gene_available:
        gt = variants[rsid].genotype   # (allele_a, allele_b)
        hap_a[rsid] = gt[0]
        hap_b[rsid] = gt[1]

    best_a = _best_match(hap_a, named_alleles, gene_available)
    best_b = _best_match(hap_b, named_alleles, gene_available)

    score_a = _compute_allele_score(
        next(na for na in named_alleles if na.name == best_a), gene_available
    )
    score_b = _compute_allele_score(
        next(na for na in named_alleles if na.name == best_b), gene_available
    )
    total_score = score_a + score_b

    # Confidence = доля покрытых defining positions
    max_possible = sum(
        _compute_allele_score(na, gene_available | set(missing_positions))
        for na in named_alleles
        if na.name in (best_a, best_b)
    )
    confidence = total_score / max_possible if max_possible > 0 else 0.0

    return DiplotypeCall(
        gene=gene,
        allele_a=best_a,
        allele_b=best_b,
        phased=True,
        confidence=round(confidence, 3),
        diplotype_score=total_score,
        missing_positions=missing_positions,
        low_quality_positions=low_qual,
    )


def _best_match(
    haplotype: dict[str, str],
    named_alleles: list[NamedAlleleDefinition],
    available_rsids: set[str],
) -> str:
    """
    Для одного хаплотипа выбирает наиболее специфичный совпавший аллель.
    Если ни один не совпал → возвращает '*1' (reference) по умолчанию.
    """
    candidates: list[tuple[int, str]] = []  # (score, name)

    for na in named_alleles:
        if _match_single_haplotype(haplotype, na, available_rsids):
            score = _compute_allele_score(na, available_rsids)
            candidates.append((score, na.name))

    if not candidates:
        logger.warning(f"Ни один аллель не совпал → fallback на *1")
        return "*1"

    # Выбираем аллель с максимальным скором (больше defining variants = специфичнее)
    candidates.sort(key=lambda x: x[0], reverse=True)
    return candidates[0][1]


def _call_unphased(
    gene: str,
    named_alleles: list[NamedAlleleDefinition],
    variants: dict[str, VariantCall],
    gene_available: set[str],
    missing_positions: list[str],
    low_qual: list[str],
) -> DiplotypeCall:
    """
    Вызов диплотипа для НЕФАЗИРОВАННЫХ данных.

    Ключевая идея: для гетерозиготных позиций мы не знаем какой аллель
    на какой хромосоме. Генерируем все возможные распределения и
    ищем пару (allele_A, allele_B) с максимальным суммарным скором,
    которая «объясняет» наблюдаемые генотипы без противоречий.
    """
    # Разбиваем варианты на гомо- и гетерозиготные
    het_rsids = [
        r for r in gene_available if variants[r].genotype[0] != variants[r].genotype[1]
    ]
    hom_rsids = [
        r for r in gene_available if variants[r].genotype[0] == variants[r].genotype[1]
    ]

    # Для гомозигот оба хаплотипа имеют одинаковый нуклеотид
    base_hap: dict[str, str] = {r: variants[r].genotype[0] for r in hom_rsids}

    # Ограничение: не более 2^12 = 4096 комбинаций для скорости
    MAX_HET = 12
    if len(het_rsids) > MAX_HET:
        logger.warning(
            f"{gene}: слишком много гетерозигот ({len(het_rsids)}), "
            f"ограничиваем до {MAX_HET} наиболее значимых"
        )
        # Оставляем только позиции, определяющие максимальное число аллелей
        importance = {
            r: sum(1 for na in named_alleles if r in na.defining_variants)
            for r in het_rsids
        }
        het_rsids = sorted(het_rsids, key=lambda r: importance[r], reverse=True)[:MAX_HET]

    best_diplotype: tuple[str, str] | None = None
    best_score = -1
    all_candidates: list[tuple[str, str, int]] = []

    # Перебираем все 2^N распределений гетерозигот по хромосомам
    for assignment in itertools.product([0, 1], repeat=len(het_rsids)):
        # assignment[i] = 0 → genotype[0] на хромосоме A
        # assignment[i] = 1 → genotype[1] на хромосоме A

        hap_a = dict(base_hap)
        hap_b = dict(base_hap)

        for rsid, side in zip(het_rsids, assignment):
            gt = variants[rsid].genotype
            hap_a[rsid] = gt[side]
            hap_b[rsid] = gt[1 - side]

        # Матчим каждый хаплотип
        matched_a: list[tuple[int, str]] = []
        matched_b: list[tuple[int, str]] = []

        for na in named_alleles:
            if _match_single_haplotype(hap_a, na, gene_available):
                matched_a.append((_compute_allele_score(na, gene_available), na.name))
            if _match_single_haplotype(hap_b, na, gene_available):
                matched_b.append((_compute_allele_score(na, gene_available), na.name))

        if not matched_a or not matched_b:
            continue

        matched_a.sort(reverse=True)
        matched_b.sort(reverse=True)

        # Проверяем consistency: диплотип должен «объяснить» все наблюдаемые варианты
        # (детальная проверка: не допускаем два одинаковых ref-аллеля объяснять разные GT)
        candidate_a = matched_a[0][1]
        candidate_b = matched_b[0][1]
        score = matched_a[0][0] + matched_b[0][0]

        diplotype_key = tuple(sorted([candidate_a, candidate_b]))
        all_candidates.append((candidate_a, candidate_b, score))

        if score > best_score:
            best_score = score
            best_diplotype = (candidate_a, candidate_b)

    if best_diplotype is None:
        logger.warning(f"{gene}: не удалось вызвать диплотип → fallback *1/*1")
        best_diplotype = ("*1", "*1")
        best_score = 0

    # Нормализуем: аллель с меньшим номером — первый (конвенция CPIC)
    a, b = best_diplotype
    if a > b:
        a, b = b, a

    # Confidence: насколько хорошо покрыты позиции
    max_possible_score = sum(
        len(na.defining_variants) for na in named_alleles if na.name in (a, b)
    )
    confidence = best_score / max_possible_score if max_possible_score > 0 else 0.0

    # Убираем дубли из all_candidates и оставляем топ-5
    seen = set()
    unique_candidates = []
    for c in sorted(all_candidates, key=lambda x: x[2], reverse=True):
        key = tuple(sorted([c[0], c[1]]))
        if key not in seen:
            seen.add(key)
            unique_candidates.append(c)
        if len(unique_candidates) == 5:
            break

    return DiplotypeCall(
        gene=gene,
        allele_a=a,
        allele_b=b,
        phased=False,
        confidence=round(confidence, 3),
        diplotype_score=best_score,
        missing_positions=missing_positions,
        low_quality_positions=low_qual,
        all_candidates=unique_candidates,
    )