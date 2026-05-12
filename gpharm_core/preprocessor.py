# gpharm_core/preprocessor.py

import subprocess
import logging
import shutil
from pathlib import Path
from dataclasses import dataclass

logger = logging.getLogger(__name__)


@dataclass
class PreprocessResult:
    vcf_path: Path          # финальный PGx-only VCF
    genome_build: str       # 'GRCh38' или 'hg19' (исходный)
    liftover_applied: bool  # был ли применён liftover
    n_variants_total: int   # вариантов после нормализации
    n_variants_pgx: int     # вариантов после фильтрации по BED


# ─────────────────────────────────────────────────────────────────────────────


def _check_tool(tool: str) -> None:
    """Проверяет что инструмент установлен, иначе сразу падает с понятной ошибкой."""
    if shutil.which(tool) is None:
        raise EnvironmentError(
            f"Инструмент '{tool}' не найден в PATH. "
            f"Установи его: conda install -c bioconda {tool}"
        )


def _run(cmd: list[str], step_name: str) -> subprocess.CompletedProcess:
    """Запускает команду, логирует её, бросает исключение при ненулевом коде возврата."""
    logger.info(f"[{step_name}] Запуск: {' '.join(cmd)}")
    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if result.returncode != 0:
        logger.error(f"[{step_name}] STDERR:\n{result.stderr}")
        raise RuntimeError(
            f"Шаг '{step_name}' завершился с ошибкой (код {result.returncode}).\n"
            f"Команда: {' '.join(cmd)}"
        )
    return result


def _count_variants(vcf_path: Path) -> int:
    """Считает число вариантов в VCF (без заголовков)."""
    result = subprocess.run(
        ["bcftools", "view", "--no-header", "-H", str(vcf_path)],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
    )
    return result.stdout.count("\n")


def _detect_genome_build(vcf_path: Path) -> str:
    """
    Пробует определить сборку из заголовка VCF.
    Ищет строки ##reference= или ##contig= с ключевыми словами.
    Возвращает 'GRCh38' или 'hg19' или 'unknown'.
    """
    hg38_markers = {"GRCh38", "hg38", "GCA_000001405.15"}
    hg19_markers = {"GRCh37", "hg19", "GCA_000001405.1", "b37", "hs37"}

    with open(vcf_path, "r") as f:
        for line in f:
            if not line.startswith("##"):
                break  # вышли из заголовка
            for marker in hg38_markers:
                if marker in line:
                    return "GRCh38"
            for marker in hg19_markers:
                if marker in line:
                    return "hg19"

    logger.warning(
        "Не удалось определить сборку из заголовка VCF. "
        "Считаем GRCh38 по умолчанию — передай build='hg19' явно если нужен liftover."
    )
    return "unknown"


def _liftover_hg19_to_hg38(
    vcf_path: Path,
    output_dir: Path,
    chain_file: Path,
    ref_hg38: Path,
) -> Path:
    """
    Запускает CrossMap для конвертации hg19 → GRCh38.

    Требования:
        pip install CrossMap
        chain_file: hg19ToHg38.over.chain.gz  (скачать с UCSC)
        ref_hg38:   GRCh38.fa (индексированный, .fai должен быть рядом)

    Возвращает путь к lifted-over VCF.
    """
    _check_tool("CrossMap.py")

    lifted_vcf = output_dir / (vcf_path.stem + ".hg38.vcf")
    unmap_file = output_dir / (vcf_path.stem + ".hg38.unmap.vcf")

    _run(
        [
            "CrossMap.py", "vcf",
            str(chain_file),
            str(vcf_path),
            str(ref_hg38),
            str(lifted_vcf),
        ],
        step_name="liftover_hg19→hg38",
    )

    # CrossMap не сообщает об unmapped вариантах через exit code —
    # логируем их число вручную.
    if unmap_file.exists():
        n_unmap = _count_variants(unmap_file)
        if n_unmap > 0:
            logger.warning(
                f"Liftover: {n_unmap} вариантов не удалось перевести в GRCh38 "
                f"(сохранены в {unmap_file})."
            )

    return lifted_vcf


def _normalize_vcf(
    vcf_path: Path,
    output_dir: Path,
    ref_hg38: Path,
) -> Path:
    """
    bcftools norm:
      -m -any   → разбивает multi-allelic сайты на отдельные строки
      -f ref    → left-aligns и нормализует indels относительно референса
      -D        → удаляет дубликаты после split
    """
    _check_tool("bcftools")

    normalized_vcf = output_dir / (vcf_path.stem + ".norm.vcf")

    _run(
        [
            "bcftools", "norm",
            "-m", "-any",       # split multi-allelic
            "-f", str(ref_hg38),  # left-align indels
            "-D",               # remove exact duplicates
            "-o", str(normalized_vcf),
            "-O", "v",          # output: plain VCF (не bgzipped)
            str(vcf_path),
        ],
        step_name="bcftools_norm",
    )

    return normalized_vcf


def _filter_pgx_regions(
    vcf_path: Path,
    output_dir: Path,
    bed_path: Path,
) -> Path:
    """
    bcftools view -R: оставляет только варианты, попадающие в PGx BED-регионы.
    Это и есть «сокращённый датасет» — всё остальное выбрасывается.
    """
    pgx_vcf = output_dir / (vcf_path.stem.replace(".norm", "") + ".pgx_only.vcf")

    _run(
        [
            "bcftools", "view",
            "-R", str(bed_path),  # фильтр по BED
            "--no-update",        # не пересчитывать AC/AN после фильтрации
            "-o", str(pgx_vcf),
            "-O", "v",
            str(vcf_path),
        ],
        step_name="bcftools_view_pgx_regions",
    )

    return pgx_vcf


# ─────────────────────────────────────────────────────────────────────────────
# Главная функция
# ─────────────────────────────────────────────────────────────────────────────


def preprocess_vcf(
    vcf_path: str | Path,
    ref_hg38: str | Path,
    bed_path: str | Path,
    output_dir: str | Path,
    chain_file: str | Path | None = None,
    force_build: str | None = None,
) -> PreprocessResult:
    """
    Полный препроцессинг VCF для PGx-пайплайна.

    Шаги:
        1. Автодетект или ручное задание сборки генома
        2. Liftover hg19 → GRCh38 (если нужно, через CrossMap)
        3. Нормализация: split multi-allelic + left-align indels (bcftools norm)
        4. Фильтрация по PGx BED-регионам (bcftools view -R)

    Аргументы:
        vcf_path    : путь к входному VCF (plain или .gz)
        ref_hg38    : путь к GRCh38.fa (должен быть индексирован .fai)
        bed_path    : путь к pgx_regions_grch38.bed
        output_dir  : куда складывать промежуточные и финальный файлы
        chain_file  : путь к hg19ToHg38.over.chain.gz (нужен только при hg19)
        force_build : 'GRCh38' или 'hg19' — принудительно задать сборку,
                      минуя автодетект (полезно если заголовок VCF неполный)

    Возвращает:
        PreprocessResult с путём к финальному pgx_only.vcf и метриками QC
    """
    vcf_path = Path(vcf_path)
    ref_hg38 = Path(ref_hg38)
    bed_path = Path(bed_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Валидация входных файлов
    for p, name in [(vcf_path, "vcf_path"), (ref_hg38, "ref_hg38"), (bed_path, "bed_path")]:
        if not p.exists():
            raise FileNotFoundError(f"Файл не найден: {name} = {p}")

    # ── Шаг 1: определение сборки ──────────────────────────────────────────
    build = force_build or _detect_genome_build(vcf_path)
    logger.info(f"Определена сборка генома: {build}")

    current_vcf = vcf_path
    liftover_applied = False

    # ── Шаг 2: liftover (только если hg19) ────────────────────────────────
    if build in ("hg19", "GRCh37"):
        if chain_file is None:
            raise ValueError(
                "Обнаружен hg19-VCF, но chain_file не передан.\n"
                "Скачай цепочку: "
                "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
            )
        chain_file = Path(chain_file)
        if not chain_file.exists():
            raise FileNotFoundError(f"chain_file не найден: {chain_file}")

        logger.info("Запускаю liftover hg19 → GRCh38...")
        current_vcf = _liftover_hg19_to_hg38(
            vcf_path=current_vcf,
            output_dir=output_dir,
            chain_file=chain_file,
            ref_hg38=ref_hg38,
        )
        liftover_applied = True

    # ── Шаг 3: нормализация ───────────────────────────────────────────────
    logger.info("Нормализую VCF (bcftools norm)...")
    normalized_vcf = _normalize_vcf(
        vcf_path=current_vcf,
        output_dir=output_dir,
        ref_hg38=ref_hg38,
    )
    n_total = _count_variants(normalized_vcf)
    logger.info(f"После нормализации: {n_total} вариантов")

    # ── Шаг 4: фильтрация по PGx-регионам ────────────────────────────────
    logger.info("Фильтрую по PGx BED-регионам...")
    pgx_vcf = _filter_pgx_regions(
        vcf_path=normalized_vcf,
        output_dir=output_dir,
        bed_path=bed_path,
    )
    n_pgx = _count_variants(pgx_vcf)
    logger.info(
        f"PGx-фильтрация: {n_total} → {n_pgx} вариантов "
        f"({100 * n_pgx / max(n_total, 1):.1f}% осталось)"
    )

    return PreprocessResult(
        vcf_path=pgx_vcf,
        genome_build=build,
        liftover_applied=liftover_applied,
        n_variants_total=n_total,
        n_variants_pgx=n_pgx,
    )
