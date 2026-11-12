import pandas as pd
import os

def calculate_gene_scores(annotated_tsv_path, positive_hpo, negative_hpo, hpo_db_path="data/genes_to_phenotype.txt"):
    # Проверяем, существует ли база HPO
    if not os.path.exists(hpo_db_path):
        print(f"ОШИБКА: Файл базы HPO не найден по пути {hpo_db_path}")
        print("Скачай его и положи в папку data/")
        return None

    print("Загружаем базу HPO... (это может занять пару секунд)")
    # Читаем базу HPO (пропускаем строки-комментарии, начинающиеся с #)
    hpo_db = pd.read_csv(hpo_db_path, sep="\t", comment="#")
    
    # Оставляем только нужные колонки: имя гена и HPO-термин
    # В базе HPO колонки называются 'gene_symbol' и 'hpo_id'
    gene_to_hpo = hpo_db.groupby('gene_symbol')['hpo_id'].apply(set).to_dict()
    print(f"База загружена. В ней {len(gene_to_hpo)} генов.")
    
    print("Загружаем аннотированные варианты...")
    df = pd.read_csv(annotated_tsv_path, sep="\t")
    
    # Оставляем только варианты, для которых мы нашли ген
    df_genes = df[df['gene'] != 'NA'].copy()
    print(f"Найдено {len(df_genes)} вариантов с генами.")
    
    # Считаем баллы для каждого гена
    results = []
    for _, row in df_genes.iterrows():
        gene = row['gene']
        
        # Если гена нет в базе HPO, даем 0 баллов
        if gene not in gene_to_hpo:
            score = 0
            matches = 0
            excluded_matches = 0
        else:
            gene_hpo_set = gene_to_hpo[gene]
            
            # 1. Плюс за совпадение с симптомами пациента
            matches = len(gene_hpo_set.intersection(positive_hpo))
            
            # 2. Минус за совпадение с ИСКЛЮЧЕННЫМИ симптомами
            excluded_matches = len(gene_hpo_set.intersection(negative_hpo))
            
            # Формула: +1 за каждый правильный симптом, -0.5 за каждый неправильный
            score = matches - (excluded_matches * 0.5)
        
        results.append({
            'gene': gene,
            'variant': f"{row['chrom']}:{row['pos']}:{row['ref']}:{row['alt']}",
            'consequence': row['consequence'],
            'score': score,
            'matches': matches,
            'excluded_matches': excluded_matches
        })
    
    # Превращаем в DataFrame и сортируем по убыванию баллов
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values(by='score', ascending=False)
    
    return results_df

if __name__ == "__main__":
    # Импортируем парсер фенотипа из нашего прошлого файла
    try:
        from phenotype_parser import parse_phenotype
    except ImportError:
        # Если запускаем не из папки src, добавляем путь
        import sys
        sys.path.append('src')
        from phenotype_parser import parse_phenotype
    
    print("--- ЗАПУСК АЛГОРИТМА СКОРИНГА ---")
    
    # 1. Загружаем симптомы
    pos_hpo, neg_hpo = parse_phenotype("data/case_01_phenotype.json")
    
    # 2. Считаем баллы
    ranked_genes = calculate_gene_scores(
        annotated_tsv_path="data/case_01_annotated.tsv",
        positive_hpo=pos_hpo,
        negative_hpo=neg_hpo
    )
    
    # 3. Показываем топ-5 результатов
    if ranked_genes is not None and not ranked_genes.empty:
        print("\n--- ТОП-5 ГЕНОВ-КАНДИДАТОВ ---")
        # Настраиваем pandas, чтобы он показывал все колонки
        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', 1000)
        print(ranked_genes.head(5).to_string(index=False))
    else:
        print("Не удалось получить результаты.")