from typing import List
import pandas as pd
from models.diplotype import DiplotypeCall


class PhenotypeEngine:
    def __init__(self, rules_star_tsv_path: str):
        """
        Инициализация напрямую из таблицы правил star-аллелей.
        """
        df = pd.read_csv(rules_star_tsv_path, sep='\t')

        # Создаем словарь быстрого перевода: (gene, diplotype) -> phenotype
        self.phenotype_map = {}
        for _, row in df.iterrows():
            if 'phenotype' in row and pd.notna(row['phenotype']):
                self.phenotype_map[(row['gene'], row['diplotype'])] = row['phenotype']

    def assign_phenotypes(self, diplotype_calls: List[DiplotypeCall]) -> List[DiplotypeCall]:
        for call in diplotype_calls:
            # Ищем фенотип в словаре
            key = (call.gene, call.diplotype)
            call.phenotype = self.phenotype_map.get(key, "Normal Metabolizer")
            call.activity_score = -1.0  # Индикатор прямого маппинга без баллов
        return diplotype_calls