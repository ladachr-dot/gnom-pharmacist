# pgx_core/hla_caller.py
"""
HLA-типирование — диспетчер для HLA-HD (WGS) и OptiType (WES/array).

Клинически значимые аллели для PGx (реакции гиперчувствительности):
    HLA-B*57:01  → абакавир       (DRESS/HSR, CPIC Level A)
    HLA-B*58:01  → аллопуринол   (SJS/TEN,   CPIC Level A)
    HLA-B*15:02  → карбамазепин  (SJS/TEN,   CPIC Level A, азиатские популяции)
    HLA-A*31:01  → карбамазепин  (DRESS,      CPIC Level A, европеоиды)

Почему два инструмента:
    HLA-HD  — поддерживает полный класс I+II, 6-digit точность, использует IPD-IMGT/HLA DB
    OptiType — только класс I (HLA-A/B/C), 4-digit, быстрый, хорошо работает на WES

Установка:
    HLA-HD:  https://www.genome.med.kyoto-u.ac.jp/HLA-HD/
             (требует bowtie2 + IPD-IMGT/HLA словарь)
    OptiType: pip install OptiType
              (требует razers3 или yara для pre-filtering)
"""

import csv
import logging
import shutil
import subprocess
import re
from enum import Enum
from pathlib import Path
from dataclasses import dataclass, field

logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# Клинически значимые аллели — единственное что нас реально интересует
# ─────────────────────────────────────────────────────────────────────────────

# Структура: {нормализованный аллель: (препарат, тип реакции, уровень доказательности)}
CLINICAL_HLA_TABLE: dict[str, tuple[str, str, str]] = {
    "HLA-B*57:01": ("Абакавир",       "HSR/DRESS",  "CPIC-A"),
    "HLA-B*58:01": ("Аллопуринол",    "SJS/TEN",    "CPIC-A"),
    "HLA-B*15:02": ("Карбамазепин",   "SJS/TEN",    "CPIC-A"),
    "HLA-A*31:01": ("Карбамазепин",   "DRESS",      "CPIC-A"),
    # Расширенный список (CPIC B / PharmGKB 1A) — менее строгие, но клинически важные
    "HLA-B*44:03": ("Ванкомицин",     "RED-man",    "PharmGKB-1A"),
    "HLA-B*35:01": ("Мирабегрон",     "DRESS",      "PharmGKB-2A"),
    "HLA-A*02:01": ("Флуклоксациллин","DILI",       "PharmGKB-2A"),
}


class HLARiskStatus(Enum):
    POSITIVE  = "POSITIVE"    # аллель обнаружен → риск реакции
    NEGATIVE  = "NEGATIVE"    # аллель не обнаружен
    UNCERTAIN = "UNCERTAIN"   # типирование не удалось / низкое качество


@dataclass
class HLARisk:
    """Клинический результат для одного аллеля."""
    allele: str            # например 'HLA-B*57:01'
    status: HLARiskStatus
    drug: str              # препарат из CLINICAL_HLA_TABLE
    reaction_type: str     # 'SJS/TEN', 'DRESS', 'HSR'
    evidence_level: str    # 'CPIC-A', 'PharmGKB-1A', ...
    note: str = ""         # дополнительное пояснение


@dataclass
class HLATypingResult:
    """Полный результат HLA-типирования образца."""
    sample: str
    method: str                                   # 'hlahd' | 'optitype'
    raw_alleles: dict[str, list[str]]             # {локус: [аллель1, аллель2]}
    clinical_risks: list[HLARisk]                 # только клинически значимые
    resolution: str                               # '4-digit' | '6-digit' | '2-digit'
    warnings: list[str] = field(default_factory=list)

    def has_risk(self, allele: str) -> bool:
        return any(r.allele == allele and r.status == HLARiskStatus.POSITIVE
                   for r in self.clinical_risks)


# ─────────────────────────────────────────────────────────────────────────────
# Утилиты
# ─────────────────────────────────────────────────────────────────────────────

def _check_tool(tool: str) -> bool:
    return shutil.which(tool) is not None


def _normalize_hla_allele(raw: str, locus_hint: str = "") -> str:
    """
    Нормализует аллель к формату HLA-X*YY:ZZ (4-digit).

    Примеры:
        'A*31:01:02:01'  → 'HLA-A*31:01'
        'B57:01'         → 'HLA-B*57:01'
        'HLA-B*57:01:01' → 'HLA-B*57:01'
        'A*02'           → 'HLA-A*02'     (2-digit, оставляем как есть)
        'Not typed'      → ''
    """
    raw = raw.strip()
    if not raw or raw.lower() in ("not typed", "na", "-", "failed", "none"):
        return ""

    # Убираем 'HLA-' префикс если есть (добавим заново)
    raw = re.sub(r"^HLA-", "", raw, flags=re.IGNORECASE)

    # Если нет '*' — добавляем (HLA-HD иногда пишет 'B57:01')
    if "*" not in raw and ":" in raw:
        locus = re.match(r"^([A-Z]+)", raw)
        number = re.sub(r"^[A-Z]+", "", raw)
        if locus:
            raw = f"{locus.group()}*{number}"

    # Добавляем локус из подсказки если он отсутствует в строке
    if locus_hint and not re.match(r"^[A-Z]", raw):
        raw = f"{locus_hint}*{raw}"

    # Обрезаем до 4-digit (первые два поля: XX:YY)
    parts = raw.split(":")
    if len(parts) >= 2:
        field1 = parts[0]   # например 'B*57'
        field2 = parts[1]   # например '01'
        normalized = f"HLA-{field1}:{field2}"
    else:
        normalized = f"HLA-{raw}"

    return normalized


def _assess_clinical_risks(
    raw_alleles: dict[str, list[str]]
) -> list[HLARisk]:
    """
    Проверяет все вызванные аллели против CLINICAL_HLA_TABLE.
    Возвращает список HLARisk для каждого клинически значимого аллеля.
    """
    risks: list[HLARisk] = []

    # Все вызванные аллели в нормализованном виде (flat list)
    called_normalized = set()
    for locus, alleles in raw_alleles.items():
        for a in alleles:
            norm = _normalize_hla_allele(a, locus_hint=locus)
            if norm:
                called_normalized.add(norm)

    for target_allele, (drug, reaction, evidence) in CLINICAL_HLA_TABLE.items():
        if target_allele in called_normalized:
            status = HLARiskStatus.POSITIVE
            note = f"Аллель {target_allele} обнаружен → риск {reaction} при применении {drug}"
            logger.warning(f"⚠️  КЛИНИЧЕСКИЙ РИСК: {note}")
        else:
            status = HLARiskStatus.NEGATIVE
            note = ""

        risks.append(HLARisk(
            allele=target_allele,
            status=status,
            drug=drug,
            reaction_type=reaction,
            evidence_level=evidence,
            note=note,
        ))

    return risks


# ─────────────────────────────────────────────────────────────────────────────
# HLA-HD (WGS — класс I + II, 6-digit)
# ─────────────────────────────────────────────────────────────────────────────

def _extract_reads_for_hla(bam_path: Path, output_dir: Path) -> tuple[Path, Path]:
    """
    Извлекает reads из MHC-региона (chr6:28,510,120-33,480,577) в FASTQ.
    Это ускоряет HLA-HD в 5-10x по сравнению с полным BAM.
    """
    r1 = output_dir / "hla_reads_R1.fastq"
    r2 = output_dir / "hla_reads_R2.fastq"

    subprocess.run([
        "samtools", "view", "-b",
        "-o", str(output_dir / "hla_region.bam"),
        str(bam_path),
        "chr6:28510120-33480577"
    ], check=True, capture_output=True)

    subprocess.run([
        "samtools", "sort", "-n",
        "-o", str(output_dir / "hla_region_sorted.bam"),
        str(output_dir / "hla_region.bam")
    ], check=True, capture_output=True)

    subprocess.run([
        "samtools", "fastq",
        "-1", str(r1), "-2", str(r2),
        "-0", "/dev/null", "-s", "/dev/null",
        str(output_dir / "hla_region_sorted.bam")
    ], check=True, capture_output=True)

    return r1, r2


def _run_hlahd(
    bam_path: Path,
    output_dir: Path,
    hlahd_dir: Path,
    hla_db_path: Path,
    threads: int,
) -> dict[str, list[str]]:
    """
    Запускает HLA-HD и парсит результаты.

    Команда:
        hlahd.sh -t {threads} -m 100 -f {freq_data} \\
            {R1.fastq} {R2.fastq} {hla_gene_split_file} {hla_db} \\
            {sample_name} {output_dir}

    Выходной файл: {output_dir}/{sample}/result/{sample}_final.result.txt
    Формат:
        HLA-A   A*02:01:01  A*32:01:01
        HLA-B   B*57:01:01  B*44:02:01
        HLA-C   C*06:02:01  C*05:01:01
        HLA-DRB1 DRB1*07:01:01  DRB1*13:01:01
        ...
    """
    sample_name = bam_path.stem

    # Шаг 1: извлечь reads из MHC-региона
    r1, r2 = _extract_reads_for_hla(bam_path, output_dir)

    hlahd_sh = hlahd_dir / "bin" / "hlahd.sh"
    freq_dir = hlahd_dir / "freq_data"
    split_file = hlahd_dir / "HLA_gene.split.txt"

    if not hlahd_sh.exists():
        raise FileNotFoundError(
            f"hlahd.sh не найден: {hlahd_sh}\n"
            "Скачай HLA-HD: https://www.genome.med.kyoto-u.ac.jp/HLA-HD/"
        )

    cmd = [
        "bash", str(hlahd_sh),
        "-t", str(threads),
        "-m", "100",          # минимальная длина ридов
        "-f", str(freq_dir),  # частоты аллелей (для приоритизации)
        str(r1), str(r2),
        str(split_file),
        str(hla_db_path),
        sample_name,
        str(output_dir),
    ]

    logger.info(f"Запуск HLA-HD: {' '.join(cmd)}")
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    if proc.returncode != 0:
        raise RuntimeError(
            f"HLA-HD завершился с ошибкой (код {proc.returncode}).\n"
            f"STDERR: {proc.stderr[:800]}"
        )

    # Парсинг результатов
    result_file = output_dir / sample_name / "result" / f"{sample_name}_final.result.txt"
    if not result_file.exists():
        raise FileNotFoundError(
            f"HLA-HD result не найден: {result_file}\n"
            f"STDOUT:\n{proc.stdout[:500]}"
        )

    return _parse_hlahd_result(result_file)


def _parse_hlahd_result(result_file: Path) -> dict[str, list[str]]:
    """
    Парсит {sample}_final.result.txt от HLA-HD.
    Каждая строка: GENE  allele1  allele2  (tab-separated)
    'Not typed' означает что аллель не удалось определить.
    """
    alleles: dict[str, list[str]] = {}
    warnings_log = []

    with open(result_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 3:
                continue

            locus = parts[0].replace("HLA-", "")  # 'A', 'B', 'DRB1', ...
            a1_raw = parts[1]
            a2_raw = parts[2]

            a1 = _normalize_hla_allele(a1_raw, locus_hint=locus)
            a2 = _normalize_hla_allele(a2_raw, locus_hint=locus)

            if not a1 and not a2:
                warnings_log.append(f"{locus}: Not typed")
                continue

            alleles[locus] = [x for x in [a1, a2] if x]

    if warnings_log:
        logger.warning(f"HLA-HD: не типированы локусы: {warnings_log}")

    return alleles


# ─────────────────────────────────────────────────────────────────────────────
# OptiType (WES — только класс I HLA-A/B/C, 4-digit)
# ─────────────────────────────────────────────────────────────────────────────

def _run_optitype(
    bam_path: Path,
    output_dir: Path,
    hla_ref_fasta: Path,
    threads: int,
) -> dict[str, list[str]]:
    """
    Запускает OptiType через pipeline:
        1. razers3 для выравнивания ридов на HLA-референс
        2. OptiType для HLA genotyping

    Выходной файл: {output_dir}/{prefix}_result.tsv
    Формат (tab-separated, первая строка — заголовок):
        \tA1\tA2\tB1\tB2\tC1\tC2\tReads\tObjective
        0\tA*02:01\tA*32:01\tB*57:01\tB*44:02\tC*06:02\tC*05:01\t1234\t0.987
    """
    sample_name = bam_path.stem
    optitype_out = output_dir / "optitype"
    optitype_out.mkdir(exist_ok=True)

    # Шаг 1: извлечь reads из HLA-региона через razers3
    hla_fastq = optitype_out / "hla_filtered.fastq"

    if not _check_tool("razers3"):
        raise EnvironmentError(
            "razers3 не найден. Установи: conda install -c bioconda razers3\n"
            "Или используй yara как альтернативу."
        )

    # Сначала BAM → FASTQ
    raw_fastq = optitype_out / "raw.fastq"
    subprocess.run(
        ["samtools", "bam2fq", "-o", str(raw_fastq), str(bam_path)],
        check=True, capture_output=True
    )

    # razers3: выравниваем только на HLA-референс (фильтрация, 95% identity)
    hla_bam = optitype_out / "hla_razers.bam"
    subprocess.run([
        "razers3", "-i", "95", "-m", "1", "-dr", "0",
        "-o", str(hla_bam),
        str(hla_ref_fasta),
        str(raw_fastq),
    ], check=True, capture_output=True)

    subprocess.run(
        ["samtools", "bam2fq", "-o", str(hla_fastq), str(hla_bam)],
        check=True, capture_output=True
    )

    # Шаг 2: OptiType
    if not _check_tool("OptiTypePipeline.py"):
        raise EnvironmentError(
            "OptiType не найден.\n"
            "Установи: pip install OptiType\n"
            "Или Docker: docker pull fred2/optitype"
        )

    cmd = [
        "OptiTypePipeline.py",
        "-i", str(hla_fastq),
        "--dna",                       # DNA-режим (не RNA)
        "-v",                          # verbose
        "-o", str(optitype_out),
        "--prefix", sample_name,
        "-t", str(threads),
    ]

    logger.info(f"Запуск OptiType: {' '.join(cmd)}")
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    if proc.returncode != 0:
        raise RuntimeError(
            f"OptiType завершился с ошибкой (код {proc.returncode}).\n"
            f"STDERR: {proc.stderr[:800]}"
        )

    return _parse_optitype_result(optitype_out, sample_name)


def _parse_optitype_result(
    optitype_out: Path,
    sample_name: str,
) -> dict[str, list[str]]:
    """
    Парсит {sample}_result.tsv от OptiType.
    Формат: строка с полями A1, A2, B1, B2, C1, C2.
    """
    # OptiType называет файл по timestamp — берём последний
    result_files = sorted(optitype_out.glob(f"*{sample_name}*result.tsv"))
    if not result_files:
        # Иногда имя файла содержит дату вместо sample_name
        result_files = sorted(optitype_out.glob("*result.tsv"))

    if not result_files:
        raise FileNotFoundError(
            f"OptiType result.tsv не найден в {optitype_out}"
        )

    result_file = result_files[-1]
    alleles: dict[str, list[str]] = {}

    with open(result_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            # Заголовки: A1, A2, B1, B2, C1, C2
            for locus in ("A", "B", "C"):
                a1_raw = row.get(f"{locus}1", "").strip()
                a2_raw = row.get(f"{locus}2", "").strip()
                a1 = _normalize_hla_allele(a1_raw, locus_hint=locus)
                a2 = _normalize_hla_allele(a2_raw, locus_hint=locus)
                if a1 or a2:
                    alleles[locus] = [x for x in [a1, a2] if x]
            break  # берём только первую строку данных

    logger.info(f"OptiType вызвал локусы: {list(alleles.keys())}")
    return alleles


# ─────────────────────────────────────────────────────────────────────────────
# Главный диспетчер
# ─────────────────────────────────────────────────────────────────────────────

def call_hla_alleles(
    bam_path: str | Path,
    output_dir: str | Path,
    input_type: str = "auto",
    hlahd_dir: str | Path | None = None,
    hla_db_path: str | Path | None = None,
    hla_ref_fasta: str | Path | None = None,
    threads: int = 4,
) -> HLATypingResult:
    """
    Диспетчер HLA-типирования: выбирает инструмент по типу данных.

    Логика выбора:
        'wgs'  → HLA-HD (класс I + II, 6-digit точность)
        'wes'  → OptiType (только класс I A/B/C, 4-digit)
        'auto' → определяет по среднему coverage MHC-региона:
                 coverage ≥ 25x  → wgs (HLA-HD)
                 coverage < 25x  → wes (OptiType)

    Аргументы:
        bam_path:       путь к BAM (GRCh38, индексирован)
        output_dir:     директория для промежуточных файлов и результатов
        input_type:     'wgs' | 'wes' | 'auto'
        hlahd_dir:      корневая папка установки HLA-HD (нужна для WGS)
        hla_db_path:    путь к IPD-IMGT/HLA словарю для HLA-HD
        hla_ref_fasta:  HLA-референс FASTA для razers3 (нужен для OptiType)
        threads:        число потоков

    Возвращает:
        HLATypingResult с raw_alleles и clinical_risks
    """
    bam_path = Path(bam_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not bam_path.exists():
        raise FileNotFoundError(f"BAM не найден: {bam_path}")

    sample_name = bam_path.stem
    warnings: list[str] = []

    # ── Автодетект типа данных ────────────────────────────────────────────────
    if input_type == "auto":
        input_type = _detect_input_type(bam_path)
        logger.info(f"Автодетект типа данных: {input_type}")

    # ── Запуск нужного инструмента ────────────────────────────────────────────
    raw_alleles: dict[str, list[str]] = {}
    method: str
    resolution: str

    if input_type == "wgs":
        if not hlahd_dir or not hla_db_path:
            raise ValueError(
                "Для WGS режима (HLA-HD) необходимо передать hlahd_dir и hla_db_path.\n"
                "Скачай HLA-HD: https://www.genome.med.kyoto-u.ac.jp/HLA-HD/"
            )
        logger.info(f"Запускаю HLA-HD для {sample_name} (WGS режим)")
        raw_alleles = _run_hlahd(
            bam_path, output_dir, Path(hlahd_dir), Path(hla_db_path), threads
        )
        method = "hlahd"
        resolution = "6-digit"

    elif input_type == "wes":
        if not hla_ref_fasta:
            raise ValueError(
                "Для WES режима (OptiType) необходимо передать hla_ref_fasta.\n"
                "Скачай: https://github.com/FRED-2/OptiType/tree/master/data/hla_reference_dna.fasta"
            )
        logger.info(f"Запускаю OptiType для {sample_name} (WES режим)")
        raw_alleles = _run_optitype(
            bam_path, output_dir, Path(hla_ref_fasta), threads
        )
        method = "optitype"
        resolution = "4-digit"

        # OptiType не типирует класс II — важное предупреждение
        warnings.append(
            "OptiType (WES): типированы только HLA-A, HLA-B, HLA-C (класс I). "
            "Класс II (DRB1, DQB1 и др.) недоступен. "
            "Для полного HLA-профиля используйте WGS + HLA-HD."
        )

    else:
        raise ValueError(f"Неизвестный input_type: '{input_type}'. Допустимые: 'wgs', 'wes', 'auto'")

    # ── Клиническая оценка рисков ─────────────────────────────────────────────
    clinical_risks = _assess_clinical_risks(raw_alleles)

    positive_risks = [r for r in clinical_risks if r.status == HLARiskStatus.POSITIVE]
    if positive_risks:
        logger.warning(
            f"ОБНАРУЖЕНЫ КЛИНИЧЕСКИ ЗНАЧИМЫЕ HLA-АЛЛЕЛИ: "
            f"{[r.allele for r in positive_risks]}"
        )

    result = HLATypingResult(
        sample=sample_name,
        method=method,
        raw_alleles=raw_alleles,
        clinical_risks=clinical_risks,
        resolution=resolution,
        warnings=warnings,
    )

    _log_summary(result)
    return result


def _detect_input_type(bam_path: Path) -> str:
    """
    Определяет тип данных по среднему coverage MHC-региона.
    ≥ 25x → WGS, < 25x → WES.
    Использует samtools coverage как быстрый способ.
    """
    try:
        proc = subprocess.run(
            [
                "samtools", "coverage",
                "-r", "chr6:28510120-33480577",
                str(bam_path),
            ],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=60
        )
        for line in proc.stdout.splitlines():
            if line.startswith("chr6"):
                cols = line.split("\t")
                mean_depth = float(cols[6]) if len(cols) > 6 else 0
                detected = "wgs" if mean_depth >= 25 else "wes"
                logger.info(f"MHC coverage = {mean_depth:.1f}x → режим '{detected}'")
                return detected
    except Exception as e:
        logger.warning(f"Не удалось определить coverage: {e} → fallback 'wgs'")

    return "wgs"


def _log_summary(result: HLATypingResult) -> None:
    logger.info("─── HLA Typing Summary ───────────────────────────")
    logger.info(f"  Образец: {result.sample}")
    logger.info(f"  Метод: {result.method.upper()}  Разрешение: {result.resolution}")
    for locus, alleles in sorted(result.raw_alleles.items()):
        logger.info(f"  HLA-{locus}: {' / '.join(alleles)}")
    logger.info("  Клинические риски:")
    for risk in result.clinical_risks:
        icon = "⚠️ " if risk.status == HLARiskStatus.POSITIVE else "✓ "
        logger.info(f"  {icon} {risk.allele} → {risk.drug} ({risk.status.value})")
    logger.info("──────────────────────────────────────────────────")