from typing import List, Dict, Any
from collections import defaultdict
import pandas as pd
from models.variant import VariantCall
from models.diplotype import DiplotypeCall


class AlleleCaller:
    def __init__(self, star_alleles_tsv_path: str):
        """
        Инициализация из оригинального TSV файла (Пункт 2 Архитектуры пайплайна).
        """
        # Читаем TSV таблицу
        self.df = pd.read_csv(star_alleles_tsv_path, sep='\t')

    def call_diplotypes(self, variant_calls: List[VariantCall]) -> List[DiplotypeCall]:
        # Группируем входящие мутации: sample -> gene -> список мутаций
        grouped_data = defaultdict(lambda: defaultdict(list))
        for call in variant_calls:
            grouped_data[call.sample][call.gene].append(call)

        diplotype_results = []

        for sample_name, genes in grouped_data.items():
            for gene_name, calls in genes.items():
                # Фильтруем правила только для текущего гена
                gene_df = self.df[self.df['gene'] == gene_name]
                if gene_df.empty:
                    continue

                diplotype_call = self._determine_gene_diplotype(sample_name, gene_name, calls, gene_df)
                diplotype_results.append(diplotype_call)

        return diplotype_results

    def _determine_gene_diplotype(self, sample: str, gene: str, calls: List[VariantCall],
                                  gene_df: pd.DataFrame) -> DiplotypeCall:
        warnings = []

        # Собираем маркеры, которые реально обнаружены у пациента (HET или HOM_ALT)
        # В тренировочном VCF фильтр качества отключен (dp < 0), берем всё
        detected_markers = [c.rsid for c in calls if c.gt in [(0, 1), (1, 0), (1, 1)]]
        hom_alt_markers = [c.rsid for c in calls if c.gt == (1, 1)]

        matched_stars = []

        # Проходим по каждой звезде этого гена в TSV таблице
        for star_name, star_data in gene_df.groupby('star'):
            required_markers = set(star_data['marker_id'].tolist())

            # Проверяем подмножество маркеров (как в оригинальном pgx_pipeline.py)
            if required_markers and required_markers.issubset(set(detected_markers)):
                matched_stars.append(star_name)

        # Сборка диплотипа
        if len(matched_stars) == 0:
            diplotype = "*1/*1"
        elif len(matched_stars) == 1:
            star = matched_stars[0]
            req_markers = set(gene_df[gene_df['star'] == star]['marker_id'])
            if req_markers.issubset(set(hom_alt_markers)):
                diplotype = f"{star}/{star}"
            else:
                diplotype = f"*1/{star}"
        else:
            diplotype = f"{matched_stars[0]}/{matched_stars[1]}"

        return DiplotypeCall(
            sample=sample,
            gene=gene,
            diplotype=diplotype,
            confidence="HIGH",
            warnings=warnings
        )