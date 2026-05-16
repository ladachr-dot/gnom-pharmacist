from dataclasses import dataclass, field
from typing import List

@dataclass
class DiplotypeCall:
    sample: str
    gene: str
    diplotype: str          # Например: "*1/*2" или "*2/*17"
    phenotype: str = "Unknown"  # Будет заполнено на Этапе 3 (Phenotype Engine)
    activity_score: float = 0.0 # Будет рассчитано на Этапе 3
    confidence: str = "LOW"    # HIGH, MODERATE, LOW (QC Layer)
    warnings: List[str] = field(default_factory=list)