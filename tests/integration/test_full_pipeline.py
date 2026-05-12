# tests/integration/test_full_pipeline.py
import pytest
from pathlib import Path
from gpharm_core import (
    preprocess_vcf, extract_pgx_variants,
    call_star_alleles, compute_coverage_qc, QCStatus
)
from tests.fixtures.NA12878_expected import EXPECTED_DIPLOTYPES, EXPECTED_HLA

NA12878_VCF = Path("tests/fixtures/NA12878.vcf")
ALLELE_DEF_DIR = Path("data/allele_definitions/")
REF_HG38 = Path("/refs/GRCh38.fa")
BED_PATH = Path("data/reference/pgx_regions_grch38.bed")

@pytest.mark.integration
@pytest.mark.skipif(
    not NA12878_VCF.exists(),
    reason="NA12878 VCF не найден — скачай тестовые данные"
)
def test_full_core_pipeline(tmp_path):
    # 1. Препроцессинг
    result = preprocess_vcf(NA12878_VCF, REF_HG38, BED_PATH, tmp_path)
    assert result.vcf_path.exists()
    assert result.n_variants_pgx > 0

    # 2. Извлечение вариантов
    variants = extract_pgx_variants(result.vcf_path)
    assert len(variants) > 100, "Слишком мало PGx вариантов"

    # 3. QC
    qc = compute_coverage_qc(result.vcf_path, ALLELE_DEF_DIR)
    assert qc.overall_status != QCStatus.FAIL, f"QC FAIL: {qc.fail_genes}"

    # 4. Star-аллели
    diplotypes = call_star_alleles(variants, ALLELE_DEF_DIR)
    assert len(diplotypes) > 0

    # 5. Проверяем правильность относительно ground truth
    for gene, (expected_a, expected_b) in EXPECTED_DIPLOTYPES.items():
        if gene not in diplotypes:
            pytest.skip(f"Ген {gene} не вызван")
        call = diplotypes[gene]
        result_pair = tuple(sorted([call.allele_a, call.allele_b]))
        expected_pair = tuple(sorted([expected_a, expected_b]))
        assert result_pair == expected_pair, (
            f"{gene}: ожидали {expected_a}/{expected_b}, "
            f"получили {call.allele_a}/{call.allele_b}"
        )