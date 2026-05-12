# tests/unit/test_named_allele_matcher.py
from gpharm_core.named_allele_matcher import (
    extract_pgx_variants, _normalize_hla_allele,
    _match_single_haplotype, NamedAlleleDefinition
)

def test_extract_variants(tmp_path):
    """Проверяем что парсер VCF корректно читает GT."""
    vcf_content = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "chr10\t94942290\trs12248560\tC\tT\t.\tPASS\t.\tGT\t0/1\n"
        "chr10\t94981296\trs4244285\tG\tA\t.\tPASS\t.\tGT\t0/0\n"
    )
    vcf_path = tmp_path / "test.vcf"
    vcf_path.write_text(vcf_content)

    variants = extract_pgx_variants(vcf_path)

    assert "rs12248560" in variants
    assert variants["rs12248560"].genotype == ("C", "T")
    assert variants["rs4244285"].genotype == ("G", "G")
    assert variants["rs12248560"].phased == False

def test_match_allele_exact():
    """Точное совпадение — аллель должен матчиться."""
    haplotype = {"rs4244285": "A", "rs12248560": "C"}
    allele = NamedAlleleDefinition(
        name="*2",
        defining_variants={"rs4244285": "A"},
        score=1,
    )
    assert _match_single_haplotype(haplotype, allele, set(haplotype.keys())) == True

def test_match_allele_mismatch():
    """Несовпадение — аллель не должен матчиться."""
    haplotype = {"rs4244285": "G"}  # ref, не *2
    allele = NamedAlleleDefinition(
        name="*2",
        defining_variants={"rs4244285": "A"},
        score=1,
    )
    assert _match_single_haplotype(haplotype, allele, {"rs4244285"}) == False

def test_match_missing_position_passes():
    """Отсутствующая позиция не должна провалить матч (ref assumed)."""
    haplotype = {}  # позиции вообще нет
    allele = NamedAlleleDefinition(
        name="*2",
        defining_variants={"rs4244285": "A"},
        score=1,
    )
    # missing → пропускаем → матч проходит
    assert _match_single_haplotype(haplotype, allele, set()) == True