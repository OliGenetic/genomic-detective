import pandas as pd

# Укажи путь к файлу
file_path = "data/case_01_variants.tsv" 

def load_variants(path):
    print(f"Читаем файл: {path}")
    # Читаем TSV файл 
    df = pd.read_csv(path, sep="\t")
    print("Первые 5 строк таблицы:")
    print(df.head())
    print(f"Всего вариантов: {len(df)}")
    return df

# Если мы запускаем этот файл напрямую, выполни функцию
if __name__ == "__main__":
    df = load_variants(file_path)