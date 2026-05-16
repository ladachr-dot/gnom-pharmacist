from dataclasses import dataclass

@dataclass
class Recommendation:
    sample: str
    drug: str
    gene: str
    phenotype: str
    effect: str
    recommendation: str
    evidence_level: str