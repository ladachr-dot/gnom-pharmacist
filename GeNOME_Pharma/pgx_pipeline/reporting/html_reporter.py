import os
from typing import List, Dict
from collections import defaultdict
from jinja2 import Environment, FileSystemLoader
from models.diplotype import DiplotypeCall
from models.recommendation import Recommendation


class HTMLReporter:
    def __init__(self, templates_dir: str):
        """
        Инициализация репортера.
        templates_dir — путь к папке, где лежит report_template.html
        """
        self.env = Environment(loader=FileSystemLoader(templates_dir))
        self.template = self.env.get_template("report_template.html")

    def generate_reports(self, diplotype_calls: List[DiplotypeCall], recommendations: List[Recommendation],
                         out_dir: str):
        """
        Группирует результаты по образцам (пациентам) и создает для каждого персональный HTML-отчет.
        """
        os.makedirs(out_dir, exist_ok=True)

        # Разделяем вызовы диплотипов и рекомендации по sample_id
        sample_diplotypes = defaultdict(list)
        for call in diplotype_calls:
            sample_diplotypes[call.sample].append(call)

        sample_recs = defaultdict(list)
        for rec in recommendations:
            sample_recs[rec.sample].append(rec)

        # Генерируем отчет для каждого уникального образца
        for sample_id in sample_diplotypes.keys():
            diplotypes = sample_diplotypes[sample_id]
            recs = sample_recs.get(sample_id, [])

            # Агрегируем общий статус уверенности (QC Layer) для образца
            # Если хотя бы один вызов LOW — общий статус LOW. Если есть MODERATE — статус MODERATE.
            confidences = [d.confidence for d in diplotypes]
            if "LOW" in confidences:
                summary_confidence = "LOW"
            elif "MODERATE" in confidences:
                summary_confidence = "MODERATE"
            else:
                summary_confidence = "HIGH"

            # Собираем все предупреждения воедино
            all_warnings = []
            for d in diplotypes:
                all_warnings.extend(d.warnings)

            # Рендерим HTML, передавая все данные в переменные шаблона
            html_content = self.template.render(
                sample_id=sample_id,
                diplotypes=diplotypes,
                recommendations=recs,
                summary_confidence=summary_confidence,
                warnings=all_warnings
            )

            # Записываем файл
            report_path = os.path.join(out_dir, f"{sample_id}_pgx_report.html")
            with open(report_path, "w", encoding="utf-8") as f:
                f.write(html_content)

            print(f"[Успех] Клинический отчет для пациента {sample_id} сохранен в: {report_path}")