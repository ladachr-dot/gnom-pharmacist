# tests/unit/test_preprocessor.py
from gpharm_core.preprocessor import _detect_genome_build
from pathlib import Path

def test_detect_build_hg38(tmp_path):
    vcf = tmp_path / "test.vcf"
    vcf.write_text("##reference=GRCh38\n#CHROM\tPOS\n")
    assert _detect_genome_build(vcf) == "GRCh38"

def test_detect_build_hg19(tmp_path):
    vcf = tmp_path / "test.vcf"
    vcf.write_text("##reference=hg19\n#CHROM\tPOS\n")
    assert _detect_genome_build(vcf) == "hg19"