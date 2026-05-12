"""
gpharm_core.
"""

__version__ = "0.1.0"
__author__ = "Dev 1"

# Публичный API:
# from gpharm_core import preprocess_vcf, call_star_alleles
from .preprocessor import preprocess_vcf, PreprocessResult
from .named_allele_matcher import extract_pgx_variants, call_star_alleles, DiplotypeCall
from .call_cyp2d6 import call_cyp2d6_with_cnv, CYP2D6CallResult
from .hla_caller import call_hla_alleles, HLATypingResult, HLARiskStatus
from .qc import compute_coverage_qc, SampleQCReport, QCStatus

__all__ = [
    # Препроцессинг
    "preprocess_vcf",
    "PreprocessResult",
    # Аллели
    "extract_pgx_variants",
    "call_star_alleles",
    "DiplotypeCall",
    # CYP2D6
    "call_cyp2d6_with_cnv",
    "CYP2D6CallResult",
    # HLA
    "call_hla_alleles",
    "HLATypingResult",
    "HLARiskStatus",
    # QC
    "compute_coverage_qc",
    "SampleQCReport",
    "QCStatus",
]