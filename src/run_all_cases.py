import os
import pandas as pd
import time
import glob

# Импортируем наши функции из других файлов
from annotate_case import annotate_variants_batch
from phenotype_parser import parse_phenotype
from scorer import calculate_gene_scores

def run_pipeline():
    # Ищем все файлы с вариантами (и в development, и в hidden, если они в подпапках)
    variant_files = sorted(glob.glob("data/**/case_*_variants.tsv", recursive=True))
    
    if not variant_files:
        print("ОШИБКА: Не найдено ни одного файла case_*_variants.tsv в папке data/")
        return

    predictions = []

    for var_file in variant_files:
        # Правильно вытаскиваем ID кейса (например, case_01)
        case_id = "_".join(os.path.basename(var_file).split('_')[:2])
        
        print(f"\n{'='*40}")
        print(f"ОБРАБОТКА КЕЙСА: {case_id}")
        print(f"{'='*40}")

        # Определяем пути к файлам (ищем в той же папке, где лежит variants.tsv)
        base_dir = os.path.dirname(var_file)
        pheno_file = os.path.join(base_dir, f"{case_id}_phenotype.json")
        annotated_file = os.path.join(base_dir, f"{case_id}_annotated.tsv")

        # 1. АННОТАЦИЯ (если файла еще нет)
        if not os.path.exists(annotated_file):
            print(f"Аннотация для {case_id} не найдена. Запускаем VEP...")
            annotate_variants_batch(var_file, annotated_file)
            time.sleep(2) # Пауза, чтобы не забанил API Ensembl
        else:
            print(f"Файл аннотации для {case_id} уже существует. Пропускаем VEP.")

        # 2. ЗАГРУЗКА СИМПТОМОВ
        if not os.path.exists(pheno_file):
            print(f"ОШИБКА: Не найден файл симптомов {pheno_file}")
            predictions.append(f"{case_id}\tno_causal_variant_in_vcf")
            continue
            
        pos_hpo, neg_hpo = parse_phenotype(pheno_file)

        # 3. СКОРИНГ
        print(f"Считаем баллы для {case_id}...")
        ranked_genes = calculate_gene_scores(annotated_file, pos_hpo, neg_hpo)

        # 4. ПРИНЯТИЕ РЕШЕНИЯ (Защита от пустых результатов)
        if ranked_genes is None or ranked_genes.empty:
            print(f"Не удалось найти кандидатов для {case_id}. Возможно, нет совпадений с HPO.")
            predictions.append(f"{case_id}\tno_causal_variant_in_vcf")
            continue

        top_candidate = ranked_genes.iloc[0]
        
        # Если лучший балл слишком низкий (например, меньше 2.0), значит нет уверенного кандидата
        if top_candidate['score'] < 2.0:
            print(f"Лучший балл ({top_candidate['score']}) слишком низкий. Пишем no_causal_variant_in_vcf.")
            predictions.append(f"{case_id}\tno_causal_variant_in_vcf")
        else:
            print(f"Лучший кандидат: {top_candidate['gene']} ({top_candidate['variant']}) с баллом {top_candidate['score']}")
            predictions.append(f"{case_id}\t{top_candidate['variant']}")

    # 5. СОХРАНЕНИЕ РЕЗУЛЬТАТОВ
    output_file = "predictions.tsv"
    with open(output_file, "w") as f:
        f.write("case_id\tprediction\n")
        for line in predictions:
            f.write(f"{line}\n")
    
    print(f"\n{'='*40}")
    print(f"ГОТОВО! Все результаты сохранены в {output_file}")
    print(f"{'='*40}")

if __name__ == "__main__":
    run_pipeline()