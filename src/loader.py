import pandas as pd
import os

# ИСПРАВЛЕНО: Функция теперь универсальна и понимает и TSV, и VCF
def load_variants(path):
    print(f"Читаем файл: {path}")
    
    if path.endswith('.vcf') or path.endswith('.vcf.gz'):
        # ДОБАВЛЕНО: Простейший парсер VCF на случай, если дадут именно его
        data = []
        with open(path, 'r') as f:
            for line in f:
                if line.startswith('#'): continue
                parts = line.strip().split('\t')
                if len(parts) >= 10:
                    # Извлекаем генотип (GT), если он есть в формате (обычно первый параметр в колонке 10)
                    fmt = parts[8].split(':')
                    sample_data = parts[9].split(':')
                    gt = sample_data[0] if 'GT' in fmt else "Unknown"
                    
                    data.append({
                        'chrom': parts[0],
                        'pos': parts[1],
                        'ref': parts[3],
                        'alt': parts[4],
                        'genotype': gt
                    })
        df = pd.DataFrame(data)
    else:
        # Читаем TSV файл 
        df = pd.read_csv(path, sep="\t")
        
    print("Первые 5 строк таблицы:")
    print(df.head())
    print(f"Всего вариантов: {len(df)}")
    return df

if __name__ == "__main__":
    file_path = "data/case_01_variants.tsv" 
    if os.path.exists(file_path):
        df = load_variants(file_path)
    else:
        print("Файл не найден, создайте тестовый файл для проверки.")
