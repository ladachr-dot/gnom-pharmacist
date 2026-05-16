import os
import gzip
from typing import Tuple, List, Dict, Any
from models.variant import VariantCall


class VCFReader:
    def __init__(self, vcf_path: str):
        """
        Инициализация кроссплатформенного парсера VCF на чистом Python.
        Работает без Си-компиляторов на Windows/Linux/macOS.
        """
        if not os.path.exists(vcf_path):
            raise FileNotFoundError(f"VCF файл не найден: {vcf_path}")

        self.vcf_path = vcf_path
        self.samples = self._extract_samples()

    def _smart_open(self):
        """Автоматически определяет, сжат файл в .gz или это обычный текст."""
        with open(self.vcf_path, 'rb') as test_f:
            magic_number = test_f.read(2)
        if magic_number == b'\x1f\x8b':
            return gzip.open(self.vcf_path, 'rt', encoding='utf-8')
        else:
            return open(self.vcf_path, 'r', encoding='utf-8')

    def _extract_samples(self) -> List[str]:
        """Быстро находит имена образцов в шапке VCF."""
        with self._smart_open() as f:
            for line in f:
                if line.startswith('#CHROM'):
                    return line.strip().split('\t')[9:]
        return []

    def fetch_targeted_variants(self, target_coords: Dict[Tuple[str, str], Dict[str, Any]]) -> List[VariantCall]:
        """
        Итерирует по VCF и извлекает целевые фармакогенетические координаты.
        Полностью сохраняет фазирование, качество покрытия GQ и глубину ридов DP.
        """
        variant_calls = []

        with self._smart_open() as f:
            for line in f:
                # Пропускаем мета-информацию
                if line.startswith('#'):
                    continue

                parts = line.strip().split('\t')
                if len(parts) < 10:
                    continue

                chrom, pos, rsid, ref, alt_field = parts[0], parts[1], parts[2], parts[3], parts[4]

                # Проверяем, есть ли текущая координата в нашей базе лекарств
                coord_key = (chrom, pos)
                coord_key_no_chr = (chrom.replace('chr', ''), pos)

                var_info = None
                if coord_key in target_coords:
                    var_info = target_coords[coord_key]
                elif coord_key_no_chr in target_coords:
                    var_info = target_coords[coord_key_no_chr]

                if var_info is None:
                    continue  # Мутация не относится к нашим лекарствам, идем дальше

                # Разбираем форматную строку, чтобы понять, где лежат GT, DP, GQ
                format_fields = parts[8].split(':')

                # Ищем индексы нужных биоинформатических полей
                gt_idx = format_fields.index('GT') if 'GT' in format_fields else -1
                dp_idx = format_fields.index('DP') if 'DP' in format_fields else -1
                gq_idx = format_fields.index('GQ') if 'GQ' in format_fields else -1
                ad_idx = format_fields.index('AD') if 'AD' in format_fields else -1

                # Проходим по всем пациентам в этой строке
                for idx, sample_name in enumerate(self.samples):
                    sample_data = parts[9 + idx].split(':')

                    # Извлекаем генотип (GT)
                    gt_str = sample_data[gt_idx] if gt_idx != -1 and gt_idx < len(sample_data) else './.'

                    # Проверяем фазирование (символ '|' означает, что фаза сохранена)
                    is_phased = '|' in gt_str

                    # Переводим строковый генотип в кортеж чисел (например, '0/1' -> (0, 1))
                    clean_gt_str = gt_str.replace('|', '/').split('/')
                    try:
                        gt = tuple(int(x) for x in clean_gt_str if x != '.')
                        if not gt:
                            gt = (0, 0)  # Фоллбэк для NO_CALL
                    except ValueError:
                        gt = (0, 0)

                    # Извлекаем качественные характеристики покрытия (QC Layer)
                    try:
                        dp = int(sample_data[dp_idx]) if dp_idx != -1 and dp_idx < len(sample_data) and sample_data[
                            dp_idx] != '.' else 0
                        gq = int(sample_data[sample_data.index('GQ')]) if 'GQ' in format_fields and sample_data[
                            format_fields.index('GQ')] != '.' else 0
                    except (ValueError, IndexError):
                        dp, gq = 0, 0

                    # Извлекаем Allele Depth (AD)
                    ad = [0]
                    if ad_idx != -1 and ad_idx < len(sample_data):
                        try:
                            ad = [int(x) for x in sample_data[ad_idx].split(',') if x != '.']
                        except ValueError:
                            ad = [0]

                    variant_calls.append(VariantCall(
                        sample=sample_name,
                        gene=var_info['gene'],
                        chrom=chrom,
                        pos=int(pos),
                        rsid=rsid if rsid != '.' else var_info.get('marker_id', '.'),
                        ref=ref,
                        alt=alt_field,
                        gt=gt,
                        phased=is_phased,
                        dp=dp,
                        gq=gq,
                        ad=ad
                    ))

        return variant_calls