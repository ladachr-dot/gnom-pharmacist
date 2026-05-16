import os
from typing import List, Dict, Any, Tuple
import pandas as pd
from models.diplotype import DiplotypeCall
from models.recommendation import Recommendation


class RecommendationEngine:
    def __init__(self, rules_dir: str):
        """
        Загружает оригинальные TSV файлы правил (Пункт 3 Архитектуры пайплайна).
        """
        self.compiled_rules: Dict[Tuple[str, str], List[Dict[str, Any]]] = {}

        # Читаем таблицы
        rules_star = pd.read_csv(os.path.join(rules_dir, 'rulesstar_alleles.tsv'), sep='\t')
        rules_snp = pd.read_csv(os.path.join(rules_dir, 'rulessnp_rules.tsv'), sep='\t')

        # Компилируем правила для Star аллелей
        for _, row in rules_star.iterrows():
            key = (row['gene'], row['diplotype'])
            if key not in self.compiled_rules:
                self.compiled_rules[key] = []

            self.compiled_rules[key].append({
                'drug': row['drug'],
                'effect': row.get('effect', 'Изменение эффективности терапии'),
                'recommendation': row.get('recommendation', 'Требуется коррекция дозы'),
                'evidence_level': 'A'  # По умолчанию для CPIC Star правил
            })

    def generate_recommendations(self, diplotype_calls: List[DiplotypeCall]) -> List[Recommendation]:
        clinical_report: List[Recommendation] = []

        for call in diplotype_calls:
            lookup_key = (call.gene, call.diplotype)

            if lookup_key in self.compiled_rules:
                for rule in self.compiled_rules[lookup_key]:
                    clinical_report.append(Recommendation(
                        sample=call.sample,
                        drug=rule["drug"],
                        gene=call.gene,
                        phenotype=call.phenotype,
                        effect=rule["effect"],
                        recommendation=rule["recommendation"],
                        evidence_level=rule["evidence_level"]
                    ))
            else:
                if call.phenotype == "Normal Metabolizer":
                    clinical_report.append(Recommendation(
                        sample=call.sample,
                        drug="All Linked Drugs",
                        gene=call.gene,
                        phenotype=call.phenotype,
                        effect="Нормальная ферментативная активность.",
                        recommendation="Применять стандартную стартовую дозу согласно инструкции.",
                        evidence_level="D"
                    ))

        return clinical_report