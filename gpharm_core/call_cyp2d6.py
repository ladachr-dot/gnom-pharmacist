# gpharm_core/cyp2d6_caller.py
"""
CYP2D6 star-allele caller — обёртка над Aldy 4.

CYP2D6 требует отдельного модуля по трём причинам:
  1. Псевдоген CYP2D7 — стандартный SNP-матчинг путается в гомологичных регионах
  2. CNV (дупликации *1x2, *2x2, делеции *5) — нужна copy-number модель
  3. Фузии CYP2D6/CYP2D7 (*68, *78, *36) — только Aldy корректно их вызывает

Aldy 4 использует integer linear programming (ILP) для точного вызова.
Мы просто запускаем его через subprocess и парсим .aldy-файл.

Установка:
    pip install aldy          # CBC solver включён по умолчанию
    samtools index file.bam   # BAM должен быть индексирован
"""

import re
import logging
import subprocess
import shutil
from pathlib import Path
from dataclasses import dataclass, field

# Импортируем контракт с Dev 2 (общая модель)
from pgx_core.named_allele_matcher import DiplotypeCall

logger = logging.getLogger(__name__)

# ─────────────────────────────────────────────────────────────────────────────
# Внутренние типы
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class AldySolution:
    """Одно решение из Aldy — может быть несколько равновероятных."""
    solution_id: int
    major_diplotype: str        # например '*1/*4+*4' (raw строка из Aldy)
    allele_a: str               # нормализованный: '*1'
    allele_b: str               # нормализованный: '*4'
    cn_configuration: str       # copy-number, например '2x*1' или '3 copies'
    confidence: int             # 0–100, из Aldy stdout
    has_duplication: bool       # есть ли '+' в диплотипе (тандемная дупликация)
    has_deletion: bool          # аллель *5 (полная делеция гена)
    has_fusion: bool            # фузия CYP2D6/CYP2D7


@dataclass
class CYP2D6CallResult:
    """
    Финальный результат CYP2D6 — расширенный DiplotypeCall с CNV-полями.
    Для передачи Dev 2 используй .to_diplotype_call()
    """
    allele_a: str
    allele_b: str
    confidence: float               # нормализовано в 0–1
    copy_number: int                # общее число копий гена
    has_duplication: bool           # танdem duplication (активность > NF)
    has_deletion: bool              # *5 = полное отсутствие гена
    has_fusion: bool                # CYP2D6/CYP2D7 химерный ген
    raw_diplotype: str              # оригинальная строка из Aldy
    all_solutions: list[AldySolution] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    def to_diplotype_call(self) -> DiplotypeCall:
        """
        Конвертирует в стандартный DiplotypeCall — контракт с Dev 2.
        CNV-специфичные поля кодируются в missing_positions как мета-флаги.
        """
        flags = []
        if self.has_duplication:
            flags.append("CNV:DUPLICATION")
        if self.has_deletion:
            flags.append("CNV:DELETION_*5")
        if self.has_fusion:
            flags.append("CNV:FUSION_CYP2D7")
        if self.copy_number != 2:
            flags.append(f"CNV:COPY_NUMBER={self.copy_number}")

        return DiplotypeCall(
            gene="CYP2D6",
            allele_a=self.allele_a,
            allele_b=self.allele_b,
            phased=True,           # Aldy всегда даёт фазированный результат
            confidence=self.confidence,
            diplotype_score=-1,    # N/A для Aldy-вызовов
            missing_positions=flags,
            low_quality_positions=self.warnings,
            all_candidates=[
                (s.allele_a, s.allele_b, s.confidence)
                for s in self.all_solutions
            ],
        )


# ─────────────────────────────────────────────────────────────────────────────
# Утилиты
# ─────────────────────────────────────────────────────────────────────────────

def _check_aldy() -> None:
    if shutil.which("aldy") is None:
        raise EnvironmentError(
            "Aldy не найден в PATH.\n"
            "Установи: pip install aldy\n"
            "Проверь установку: aldy test"
        )


def _detect_sequencing_profile(bam_path: Path) -> str:
    """
    Heuristic для автоопределения профиля секвенирования по имени файла.
    Пользователь может переопределить через аргумент profile=.

    Aldy profiles:
        wgs / illumina  — WGS ≥ 40x (рекомендовано)
        wes / exome     — WES (только 2 копии, без CNV!)
        pgx1/pgx2/pgx3  — PGRNseq capture panels
        10x             — 10X Genomics linked reads
    """
    name = bam_path.name.lower()
    if any(k in name for k in ("wgs", "genome", "whole_genome")):
        return "wgs"
    if any(k in name for k in ("wes", "exome", "wxs")):
        return "wes"
    if "pgx" in name or "pgrnseq" in name:
        return "pgx2"
    # По умолчанию — WGS (наиболее распространённый в клинике)
    return "wgs"


def _normalize_allele(raw: str) -> str:
    """
    Нормализует аллель из Aldy в стандартный CPIC-формат.

    Примеры:
        '4.021'  → '*4'    (minor sub-allele → major)
        '1.016'  → '*1'
        '139'    → '*139'
        '*4+*4'  → '*4'    (tandem — берём первый мажорный)
        '78'     → '*78'

    Aldy может возвращать minor аллели вида '4.021' — нас интересует
    только major (до точки) для совместимости с CPIC phenotype tables.
    """
    raw = raw.strip().lstrip("*")
    # Убираем minor sub-allele суффикс (после точки)
    major = raw.split(".")[0]
    # Убираем ALDY-специфичные суффиксы типа '.ALDY_2'
    major = re.sub(r"\.ALDY.*$", "", major)
    # Убираем дополнительные варианты (+rs...)
    major = re.sub(r"\s*\+rs\d+.*$", "", major)
    major = re.sub(r"\s*-rs\d+.*$", "", major)
    return f"*{major}"


def _parse_aldy_output(aldy_file: Path) -> list[AldySolution]:
    """
    Парсит .aldy файл — детальный вывод Aldy с колонками:
    #Sample Gene SolutionID Major Minor Copy Allele Location Type Coverage Effect dbSNP Code Status

    Нас интересует только колонка Major (индекс 3) по первой строке каждого SolutionID.

    Пример строки:
        NA10860  CYP2D6  1  *1/*4+*4.021  1.001;4;4.021  0  1.001  ...
    """
    solutions: dict[int, AldySolution] = {}

    if not aldy_file.exists():
        raise FileNotFoundError(f"Aldy output не найден: {aldy_file}")

    with open(aldy_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith("#") or not line:
                continue

            cols = line.split("\t")
            if len(cols) < 6:
                continue

            try:
                solution_id = int(cols[2])
            except ValueError:
                continue

            if solution_id in solutions:
                continue  # уже добавили эту решение (берём только первую строку)

            major_raw = cols[3]  # например '*1/*4+*4.021'

            # Разбиваем диплотип по '/'
            parts = major_raw.lstrip("*").split("/")
            if len(parts) != 2:
                logger.warning(f"Неожиданный формат диплотипа: {major_raw}")
                continue

            allele_a_raw = parts[0]
            allele_b_raw = parts[1]

            # Определяем флаги CNV
            full_str = major_raw
            has_duplication = "+" in full_str  # тандемная дупликация: *4+*4
            has_deletion = "*5" in full_str    # *5 = полная делеция гена
            has_fusion = any(
                x in full_str for x in ("*68", "*78", "*13")
            )  # известные CYP2D6/CYP2D7 фузии

            solutions[solution_id] = AldySolution(
                solution_id=solution_id,
                major_diplotype=major_raw,
                allele_a=_normalize_allele(allele_a_raw.split("+")[0]),
                allele_b=_normalize_allele(allele_b_raw.split("+")[0]),
                cn_configuration="",  # заполняется из stdout
                confidence=100,       # заполняется из stdout
                has_duplication=has_duplication,
                has_deletion=has_deletion,
                has_fusion=has_fusion,
            )

    return list(solutions.values())


def _parse_aldy_stdout(stdout: str, solutions: list[AldySolution]) -> list[AldySolution]:
    """
    Извлекает confidence и CN-конфигурацию из stdout Aldy.

    Пример stdout:
        Potential CYP2D6 gene structures for NA07048:
          1: 2x*1 (confidence: 100%)
        Best CYP2D6 star-alleles for NA07048:
          1: *1 / *4.021 (confidence=100%)
    """
    # Confidence для каждого solution
    conf_pattern = re.compile(r"(\d+):\s+.+\(confidence[=:]?\s*(\d+)%\)")
    cn_pattern = re.compile(r"(\d+):\s+([\dx*+]+)\s+\(confidence")

    conf_map: dict[int, int] = {}
    cn_map: dict[int, str] = {}

    for line in stdout.splitlines():
        m = conf_pattern.search(line)
        if m:
            conf_map[int(m.group(1))] = int(m.group(2))

        m2 = cn_pattern.search(line)
        if m2 and "gene structure" in stdout[max(0, stdout.find(line) - 200):stdout.find(line)]:
            cn_map[int(m2.group(1))] = m2.group(2)

    for s in solutions:
        if s.solution_id in conf_map:
            s.confidence = conf_map[s.solution_id]
        if s.solution_id in cn_map:
            s.cn_configuration = cn_map[s.solution_id]

    return solutions


# ─────────────────────────────────────────────────────────────────────────────
# Главная функция
# ─────────────────────────────────────────────────────────────────────────────

def call_cyp2d6_with_cnv(
    bam_path: str | Path,
    output_dir: str | Path,
    profile: str | None = None,
    genome: str = "hg38",
    cn_neutral_region: str | None = None,
    solver: str = "any",
    cn_override: str | None = None,
) -> CYP2D6CallResult:
    """
    Запускает Aldy 4 на BAM-файле и возвращает CYP2D6 диплотип с CNV-информацией.

    Алгоритм Aldy (ILP):
        Шаг 1 — Copy Number Detection:
            Сравнивает coverage CYP2D6-региона с CN-нейтральным регионом (CYP2D8).
            Определяет число копий гена (1, 2, 3+ или делеция).

        Шаг 2 — Major Star-Allele Calling:
            Для каждой конфигурации CN перебирает возможные комбинации major аллелей.
            Использует ILP для максимизации likelihood варианта.

        Шаг 3 — Minor Star-Allele Calling (Phasing):
            Определяет точный minor аллель (*4.021 vs *4.006) и фазирует варианты.
            Для long reads использует read-cloud информацию.

    Аргументы:
        bam_path:           путь к BAM (должен быть индексирован .bai)
        output_dir:         куда писать .aldy файл
        profile:            профиль секвенирования: 'wgs'|'wes'|'pgx1'|'pgx2'|'pgx3'|'10x'
                            None → автодетект по имени файла
        genome:             'hg38' (GRCh38) или 'hg19' (GRCh37)
        cn_neutral_region:  кастомный CN-нейтральный регион 'chr:start-end'
                            None → Aldy использует CYP2D8 по умолчанию
        solver:             'any'|'cbc'|'gurobi'
        cn_override:        принудительно задать CN конфигурацию, например '2'
                            (полезно для WES где CN-detection невозможен)

    Возвращает:
        CYP2D6CallResult с лучшим решением + все альтернативы
    """
    bam_path = Path(bam_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Проверки
    _check_aldy()
    if not bam_path.exists():
        raise FileNotFoundError(f"BAM не найден: {bam_path}")

    bai_path = bam_path.with_suffix(".bam.bai")
    bai_alt = Path(str(bam_path) + ".bai")
    if not bai_path.exists() and not bai_alt.exists():
        logger.warning(
            f"BAM index (.bai) не найден для {bam_path.name}. "
            f"Создаю: samtools index {bam_path}"
        )
        subprocess.run(["samtools", "index", str(bam_path)], check=True)

    # Автодетект профиля
    seq_profile = profile or _detect_sequencing_profile(bam_path)
    logger.info(f"CYP2D6: профиль секвенирования = '{seq_profile}'")

    if seq_profile in ("wes", "exome", "wxs"):
        logger.warning(
            "WES-профиль: Aldy предполагает ровно 2 копии CYP2D6. "
            "CNV (дупликации, *5) не будут определены корректно."
        )

    # Имя выходного файла: {sample}.CYP2D6.aldy
    sample_name = bam_path.stem
    aldy_output = output_dir / f"{sample_name}.CYP2D6.aldy"

    # ── Сборка команды Aldy ───────────────────────────────────────────────────
    cmd = [
        "aldy", "genotype",
        "-g", "CYP2D6",
        "-p", seq_profile,
        "--genome", genome,
        "-o", str(aldy_output),
        "--solver", solver,
        "--multiple-warn-level", "2",   # предупреждать при нескольких optimal решениях
    ]

    if cn_neutral_region:
        cmd += ["-n", cn_neutral_region]

    if cn_override:
        cmd += ["-c", cn_override]
        logger.info(f"CYP2D6: CN принудительно задан = {cn_override}")

    cmd.append(str(bam_path))

    # ── Запуск ────────────────────────────────────────────────────────────────
    logger.info(f"Запуск Aldy: {' '.join(cmd)}")
    proc = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    # Aldy пишет результаты в stdout (не stderr), exit code 0 = успех
    if proc.returncode != 0:
        logger.error(f"Aldy STDERR:\n{proc.stderr}")
        raise RuntimeError(
            f"Aldy завершился с ошибкой (код {proc.returncode}).\n"
            f"Команда: {' '.join(cmd)}\n"
            f"STDERR: {proc.stderr[:500]}"
        )

    logger.debug(f"Aldy STDOUT:\n{proc.stdout}")

    # ── Парсинг результатов ───────────────────────────────────────────────────
    solutions = _parse_aldy_output(aldy_output)

    if not solutions:
        raise ValueError(
            f"Aldy не вернул ни одного решения для {bam_path.name}. "
            f"Проверь логи. Попробуй запустить вручную:\n{' '.join(cmd)}"
        )

    solutions = _parse_aldy_stdout(proc.stdout, solutions)

    # Лучшее решение — первое (Aldy сортирует по confidence)
    best = solutions[0]

    # Подсчёт copy number из конфигурации (например '2x*1' → 2)
    cn_match = re.search(r"(\d+)x", best.cn_configuration)
    copy_number = int(cn_match.group(1)) if cn_match else 2

    warnings: list[str] = []

    # Предупреждение о нескольких равновероятных решениях
    if len(solutions) > 1:
        alts = " | ".join(s.major_diplotype for s in solutions[1:])
        msg = f"Несколько равновероятных диплотипов: [{alts}] → требует ручной проверки"
        logger.warning(msg)
        warnings.append(msg)

    # Предупреждение о структурных вариантах
    if best.has_fusion:
        msg = f"Обнаружена фузия CYP2D6/CYP2D7 ({best.major_diplotype}) — проверь вручную"
        logger.warning(msg)
        warnings.append(msg)

    if best.has_duplication:
        logger.info(f"CYP2D6 тандемная дупликация: {best.major_diplotype}")

    if best.has_deletion:
        logger.info(f"CYP2D6 *5 делеция (полное отсутствие гена): {best.major_diplotype}")

    result = CYP2D6CallResult(
        allele_a=best.allele_a,
        allele_b=best.allele_b,
        confidence=best.confidence / 100.0,
        copy_number=copy_number,
        has_duplication=best.has_duplication,
        has_deletion=best.has_deletion,
        has_fusion=best.has_fusion,
        raw_diplotype=best.major_diplotype,
        all_solutions=solutions,
        warnings=warnings,
    )

    logger.info(
        f"CYP2D6: {result.allele_a}/{result.allele_b}  "
        f"CN={result.copy_number}  conf={result.confidence:.0%}  "
        f"dup={result.has_duplication}  del={result.has_deletion}  fusion={result.has_fusion}"
    )

    return result