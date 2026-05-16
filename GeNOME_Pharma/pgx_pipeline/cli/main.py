import os
import argparse
import json  # Используем JSON или YAML для хранения конфигураций (п. 2 и 5 ТЗ)
import pandas as pd
import sys

# Добавляем корневую директорию проекта в пути поиска модулей,
# чтобы импорты внутри папок работали корректно
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Импортируем наши новые спроектированные компоненты
from pgx_pipeline.io.vcf_reader import VCFReader
from pgx_pipeline.calling.allele_caller import AlleleCaller
from pgx_pipeline.calling.phenotype_engine import PhenotypeEngine
from pgx_pipeline.rules.recommendation_engine import RecommendationEngine
from pgx_pipeline.reporting.html_reporter import HTMLReporter


def load_target_variants(variants_tsv_path: str) -> dict:
    """
    Вспомогательный метод для загрузки целевых фармакогенетических координат
    из базы pgx_variants.tsv и приведения их к формату поиска для cyvcf2.
    """
    if not os.path.exists(variants_tsv_path):
        raise FileNotFoundError(f"Файл целевых координат не найден: {variants_tsv_path}")

    df = pd.read_csv(variants_tsv_path, sep='\t')
    target_coords = {}
    for _, row in df.iterrows():
        chrom = str(row['chrom'])
        pos = str(row['pos'])

        # Индексируем по кортежу (хромосома, позиция)
        info = {
            'gene': row['gene'],
            'marker_id': row['marker_id'],
            'ref': row['ref'],
            'alt': row['alt'],
            'marker_type': row.get('marker_type', 'SNP')
        }
        target_coords[(chrom, pos)] = info
        # Дублируем ключ без префикса 'chr' на случай разных стандартов в VCF
        target_coords[(chrom.replace('chr', ''), pos)] = info

    return target_coords


import os
import argparse
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pgx_pipeline.io.vcf_reader import VCFReader
from pgx_pipeline.calling.allele_caller import AlleleCaller
from pgx_pipeline.calling.phenotype_engine import PhenotypeEngine
from pgx_pipeline.rules.recommendation_engine import RecommendationEngine
from pgx_pipeline.reporting.html_reporter import HTMLReporter
from pgx_pipeline.cli.main import load_target_variants # Оставляем старый загрузчик pgx_variants.tsv

def main():
    parser = argparse.ArgumentParser(description="TSV-Driven PGx Pipeline")
    parser.add_argument("--vcf", required=True, help="Путь к VCF")
    parser.add_argument("--config_dir", required=True, help="Папка со всеми TSV файлами")
    parser.add_argument("--out_dir", default="./output", help="Выходная папка")
    args = parser.parse_args()

    # Нормализуем пути Windows
    config_dir = os.path.abspath(os.path.normpath(args.config_dir))
    reports_out_dir = os.path.join(os.path.abspath(os.path.normpath(args.out_dir)), "reports")

    # Пути к твоим оригинальным TSV файлам
    variants_tsv = os.path.join(config_dir, "pgx_variants.tsv")
    star_alleles_tsv = os.path.join(config_dir, "star_alleles.tsv")
    rules_star_tsv = os.path.join(config_dir, "rulesstar_alleles.tsv")
    templates_dir = os.path.join(config_dir, "templates")

    print("[1/5] Экстракция маркеров...")
    target_coords = load_target_variants(variants_tsv)
    vcf_reader = VCFReader(args.vcf)
    variant_calls = vcf_reader.fetch_targeted_variants(target_coords)

    print("[2/5] Вызов звездных аллелей из star_alleles.tsv...")
    allele_caller = AlleleCaller(star_alleles_tsv)
    diplotype_calls = allele_caller.call_diplotypes(variant_calls)

    print("[3/5] Определение фенотипов из rulesstar_alleles.tsv...")
    phenotype_engine = PhenotypeEngine(rules_star_tsv)
    diplotype_calls = phenotype_engine.assign_phenotypes(diplotype_calls)

    print("[4/5] Подбор рекомендаций из таблиц правил...")
    recommendation_engine = RecommendationEngine(config_dir)
    recommendations = recommendation_engine.generate_recommendations(diplotype_calls)

    print("[5/5] Генерация HTML через Jinja2...")
    reporter = HTMLReporter(templates_dir)
    reporter.generate_reports(diplotype_calls, recommendations, reports_out_dir)
    print("Пайплайн успешно завершил работу!")

if __name__ == "__main__":
    main()