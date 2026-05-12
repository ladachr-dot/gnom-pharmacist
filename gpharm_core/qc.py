# pgx_core/qc.py
"""
QC-модуль: оценка покрытия PGx-позиций в VCF.

Задача: для каждого гена из CPIC allele definition tables
посчитать какой % defining positions присутствует в VCF пациента.

Пороги (основаны на практике PharmCAT и CPIC):
    ≥ 90%  → PASS   — надёжный вызов
    70–89% → WARN   — вызов возможен, но с пониженным confidence
    < 70%  → FAIL   — no-call, результат ненадёжен

Дополнительно:
    - Проверяем FILTER поле: PASS vs LOW_QUAL/LOWQ/etc.
    - Проверяем критичные позиции (defining для >50% аллелей) — всегда WARN если missing
    - Генерируем общий QC-отчёт по образцу
"""

import json
import logging
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path

logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# Типы
# ─────────────────────────────────────────────────────────────────────────────

class QCStatus(Enum):
    PASS = "PASS"   # ≥ 90% позиций покрыто
    WARN = "WARN"   # 70–89% позиций покрыто
    FAIL = "FAIL"   # < 70% позиций покрыто → no-call рекомендован


@dataclass
class PositionQC:
    """QC одной позиции."""
    rsid: str
    present: bool           # есть ли в VCF
    filter_pass: bool       # FILTER == PASS / '.' (если present)
    is_critical: bool       # defining для >50% аллелей этого гена
    genotype: str = ""      # GT если присутствует


@dataclass
class GeneQCResult:
    """QC одного гена."""
    gene: str
    status: QCStatus

    n_required: int         # всего позиций в allele definition
    n_present: int          # присутствуют в VCF (любой FILTER)
    n_pass: int             # присутствуют с FILTER=PASS
    coverage_pct: float     # n_present / n_required * 100

    missing_positions: list[str]         # rsID которых нет в VCF
    low_quality_positions: list[str]     # rsID с FILTER != PASS
    critical_missing: list[str]          # критичные отсутствующие позиции

    positions: list[PositionQC] = field(default_factory=list)  # детали по каждой позиции

    @property
    def is_callable(self) -> bool:
        """Можно ли делать вызов для этого гена."""
        return self.status in (QCStatus.PASS, QCStatus.WARN)

    def summary_line(self) -> str:
        crit = f" | CRITICAL MISSING: {self.critical_missing}" if self.critical_missing else ""
        return (
            f"{self.gene}: {self.status.value} "
            f"({self.n_present}/{self.n_required} = {self.coverage_pct:.1f}%)"
            f"{crit}"
        )


@dataclass
class SampleQCReport:
    """Общий QC-отчёт по образцу."""
    sample_name: str
    vcf_path: str
    total_variants_in_vcf: int

    gene_results: dict[str, GeneQCResult]  # {gene: GeneQCResult}

    @property
    def pass_genes(self) -> list[str]:
        return [g for g, r in self.gene_results.items() if r.status == QCStatus.PASS]

    @property
    def warn_genes(self) -> list[str]:
        return [g for g, r in self.gene_results.items() if r.status == QCStatus.WARN]

    @property
    def fail_genes(self) -> list[str]:
        return [g for g, r in self.gene_results.items() if r.status == QCStatus.FAIL]

    @property
    def overall_status(self) -> QCStatus:
        """
        Общий статус образца:
            FAIL  — есть хотя бы один FAIL по клинически критичному гену
            WARN  — есть хотя бы один WARN
            PASS  — все гены PASS
        """
        # Гены с CPIC Level A — критичные, их FAIL блокирует отчёт
        critical_genes = {
            "CYP2D6", "CYP2C19", "CYP2C9", "DPYD", "TPMT", "NUDT15",
            "SLCO1B1", "VKORC1", "G6PD",
        }
        for gene in self.fail_genes:
            if gene in critical_genes:
                return QCStatus.FAIL
        if self.warn_genes:
            return QCStatus.WARN
        return QCStatus.PASS

    def to_tsv(self) -> str:
        """Генерирует TSV-строки для записи в файл."""
        lines = ["Gene\tStatus\tN_Required\tN_Present\tCoverage_pct\tMissing\tCritical_Missing"]
        for gene, r in sorted(self.gene_results.items()):
            lines.append(
                f"{gene}\t{r.status.value}\t{r.n_required}\t{r.n_present}\t"
                f"{r.coverage_pct:.1f}\t{','.join(r.missing_positions)}\t"
                f"{','.join(r.critical_missing)}"
            )
        return "\n".join(lines)


# ─────────────────────────────────────────────────────────────────────────────
# Вспомогательные функции
# ─────────────────────────────────────────────────────────────────────────────

# Значения FILTER, которые считаем «прошедшими» QC
_PASS_FILTERS = {"PASS", ".", ""}


def _load_required_positions(json_path: Path) -> tuple[str, dict[str, int]]:
    """
    Загружает из CPIC JSON список всех rsID, нужных для гена,
    и для каждого считает сколько аллелей он определяет.

    Возвращает: (gene_name, {rsid: n_alleles_it_defines})
    """
    with open(json_path) as f:
        data = json.load(f)

    gene = data["gene"]

    # rsID в том же порядке что alleles[] в namedAlleles
    variant_rsids: list[str] = [
        v.get("rsid") or v.get("chromosomeHgvsName", f"pos_{i}")
        for i, v in enumerate(data["variants"])
    ]

    # Считаем сколько named alleles задействует каждую позицию
    rsid_usage: dict[str, int] = {rsid: 0 for rsid in variant_rsids}
    for na in data["namedAlleles"]:
        for rsid, base in zip(variant_rsids, na["alleles"]):
            if base is not None:  # null = позиция не задействована
                rsid_usage[rsid] += 1

    # Убираем позиции с нулевым использованием (технические, не влияют на вызов)
    rsid_usage = {r: c for r, c in rsid_usage.items() if c > 0}

    return gene, rsid_usage


def _parse_vcf_index(vcf_path: Path) -> tuple[dict[str, tuple[str, str]], int]:
    """
    Быстро читает VCF и строит индекс {rsid: (genotype_str, filter)}.
    Возвращает также общее число вариантов.

    Не использует pysam — работает на plain VCF после preprocessor.
    """
    index: dict[str, tuple[str, str]] = {}
    total = 0

    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            total += 1
            cols = line.split("\t", 10)
            if len(cols) < 10:
                continue

            rsid   = cols[2]
            filt   = cols[6].strip()
            fmt    = cols[8]
            sample = cols[9]

            if not rsid.startswith("rs"):
                continue

            # Извлекаем GT
            fmt_fields = fmt.split(":")
            try:
                gt_idx = fmt_fields.index("GT")
            except ValueError:
                continue

            sample_fields = sample.split(":")
            if gt_idx >= len(sample_fields):
                continue

            gt_raw = sample_fields[gt_idx].strip()
            index[rsid] = (gt_raw, filt)

    return index, total


def _classify_status(coverage_pct: float) -> QCStatus:
    if coverage_pct >= 90.0:
        return QCStatus.PASS
    elif coverage_pct >= 70.0:
        return QCStatus.WARN
    else:
        return QCStatus.FAIL


# ─────────────────────────────────────────────────────────────────────────────
# Главная функция
# ─────────────────────────────────────────────────────────────────────────────

def compute_coverage_qc(
    vcf_path: str | Path,
    allele_def_dir: str | Path,
    gene: str | None = None,
    warn_threshold: float = 90.0,
    fail_threshold: float = 70.0,
    output_tsv: str | Path | None = None,
) -> SampleQCReport:
    """
    Для каждого PGx-гена считает % defining positions, присутствующих в VCF.

    Алгоритм:
        1. Загружает CPIC allele definition JSON → список required rsID
        2. Читает VCF → строит индекс {rsid: (gt, filter)}
        3. Для каждого rsID: present? filter_pass? is_critical?
        4. Считает coverage_pct = n_present / n_required * 100
        5. Классифицирует: PASS / WARN / FAIL
        6. Критичные позиции (defining для >50% аллелей) → WARN даже при общем PASS

    Аргументы:
        vcf_path:        путь к PGx-only VCF (после preprocess_vcf)
        allele_def_dir:  папка с CPIC JSON
        gene:            если задан — только этот ген; иначе все JSON в папке
        warn_threshold:  порог для WARN (default 90%)
        fail_threshold:  порог для FAIL (default 70%)
        output_tsv:      если задан — сохраняет TSV-отчёт по этому пути

    Возвращает:
        SampleQCReport
    """
    vcf_path = Path(vcf_path)
    allele_def_dir = Path(allele_def_dir)

    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF не найден: {vcf_path}")

    sample_name = vcf_path.stem.replace(".pgx_only", "").replace(".norm", "")

    # Читаем VCF один раз — строим индекс
    logger.info(f"Читаю VCF для QC: {vcf_path.name}")
    vcf_index, total_variants = _parse_vcf_index(vcf_path)
    logger.info(f"VCF содержит {total_variants} вариантов, {len(vcf_index)} с rsID")

    # Определяем какие JSON обрабатывать
    if gene:
        json_files = [allele_def_dir / f"{gene}.json"]
    else:
        json_files = sorted(allele_def_dir.glob("*.json"))

    gene_results: dict[str, GeneQCResult] = {}

    for json_path in json_files:
        if not json_path.exists():
            logger.warning(f"Файл определений не найден: {json_path}")
            continue

        gene_name, rsid_usage = _load_required_positions(json_path)
        if not rsid_usage:
            continue

        n_required = len(rsid_usage)
        n_alleles_total = max(rsid_usage.values()) if rsid_usage else 1

        # Порог «критичности»: позиция defining для >50% аллелей = критичная
        critical_threshold = n_alleles_total * 0.5

        positions: list[PositionQC] = []
        n_present = 0
        n_pass = 0
        missing: list[str] = []
        low_qual: list[str] = []
        critical_missing: list[str] = []

        for rsid, n_alleles_using in rsid_usage.items():
            is_critical = n_alleles_using >= critical_threshold

            if rsid in vcf_index:
                gt_raw, filt = vcf_index[rsid]

                # Пропускаем missing genotype (./.)
                if "." in gt_raw:
                    positions.append(PositionQC(
                        rsid=rsid, present=False,
                        filter_pass=False, is_critical=is_critical
                    ))
                    missing.append(rsid)
                    if is_critical:
                        critical_missing.append(rsid)
                    continue

                filter_pass = filt in _PASS_FILTERS
                n_present += 1
                if filter_pass:
                    n_pass += 1
                else:
                    low_qual.append(rsid)

                positions.append(PositionQC(
                    rsid=rsid, present=True,
                    filter_pass=filter_pass,
                    is_critical=is_critical,
                    genotype=gt_raw,
                ))
            else:
                # Позиции нет в VCF вообще
                positions.append(PositionQC(
                    rsid=rsid, present=False,
                    filter_pass=False, is_critical=is_critical
                ))
                missing.append(rsid)
                if is_critical:
                    critical_missing.append(rsid)

        coverage_pct = (n_present / n_required * 100) if n_required > 0 else 0.0
        status = _classify_status(coverage_pct)

        # Повышаем до WARN если есть критичные missing — даже при PASS по coverage
        if status == QCStatus.PASS and critical_missing:
            status = QCStatus.WARN
            logger.warning(
                f"{gene_name}: общее покрытие PASS ({coverage_pct:.1f}%), "
                f"но отсутствуют критичные позиции: {critical_missing} → WARN"
            )

        result = GeneQCResult(
            gene=gene_name,
            status=status,
            n_required=n_required,
            n_present=n_present,
            n_pass=n_pass,
            coverage_pct=coverage_pct,
            missing_positions=missing,
            low_quality_positions=low_qual,
            critical_missing=critical_missing,
            positions=positions,
        )

        gene_results[gene_name] = result
        logger.info(result.summary_line())

    report = SampleQCReport(
        sample_name=sample_name,
        vcf_path=str(vcf_path),
        total_variants_in_vcf=total_variants,
        gene_results=gene_results,
    )

    # Итоговый лог
    logger.info(
        f"QC Summary [{sample_name}]: "
        f"PASS={len(report.pass_genes)} "
        f"WARN={len(report.warn_genes)} "
        f"FAIL={len(report.fail_genes)} "
        f"| Overall: {report.overall_status.value}"
    )
    if report.fail_genes:
        logger.error(f"FAIL гены (no-call рекомендован): {report.fail_genes}")
    if report.warn_genes:
        logger.warning(f"WARN гены (пониженный confidence): {report.warn_genes}")

    # Сохраняем TSV если нужно
    if output_tsv:
        output_tsv = Path(output_tsv)
        output_tsv.parent.mkdir(parents=True, exist_ok=True)
        output_tsv.write_text(report.to_tsv())
        logger.info(f"QC TSV сохранён: {output_tsv}")

    return report