from dataclasses import dataclass
from typing import Tuple, List, Dict, Any


# ==========================================
# 3. NEW DATA MODEL (dataclasses)
# ==========================================

@dataclass
class VariantCall:
    sample: str
    gene: str  # Добавили связь с геном прямо в модель вызова
    chrom: str
    pos: int
    rsid: str
    ref: str
    alt: str
    gt: Tuple[int, ...]  # Теперь это кортеж чисел, например: (0, 1) для 0/1 или (1, 1) для 1/1
    phased: bool  # Флаг фазирования (True если '|', False если '/')
    dp: int  # Глубина чтения (Read Depth) в этой точке для образца
    gq: int  # Качество генотипа (Genotype Quality)
    ad: List[int]  # Аллельная глубина (Allele Depth) - покрытие референса и альтернативы
