import os
import pandas as pd
import time
import glob

from annotate_case import annotate_variants_batch
from phenotype_parser import parse_phenotype
from scorer import calculate_gene_scores  
from loader import load_variants

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

        # 1. АННОТАЦИЯ
        if not os.path.exists(annotated_file):
            print(f"Аннотация для {case_id} не найдена. Запускаем VEP...")
            # Если подали VCF, конвертируем его в TSV для аннотатора
            if var_file.endswith(('.vcf', '.vcf.gz')):
                df_vcf = load_variants(var_file)
                temp_tsv = os.path.join(base_dir, f"{case_id}_temp.tsv")
                df_vcf.to_csv(temp_tsv, sep="\t", index=False)
                annotate_variants_batch(temp_tsv, annotated_file)
                os.remove(temp_tsv)
            else:
                annotate_variants_batch(var_file, annotated_file)
            time.sleep(2) 
        else:
            print(f"Файл аннотации для {case_id} уже существует. Пропускаем VEP.")

        # 2. ЗАГРУЗКА СИМПТОМОВ
        if not os.path.exists(pheno_file):
            print(f"ОШИБКА: Не найден файл симптомов {pheno_file}")
            predictions.append(f"{case_id}\tno_causal_variant_in_vcf")
            continue
            
        pos_hpo, neg_hpo, _, _ = parse_phenotype(pheno_file)

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
            print(f"Лучший балл ({top_candidate['score']}) слишком низкий.")
            predictions.append(f"{case_id}\tno_causal_variant_in_vcf\tБаллы ниже порога достоверности")
        
        # ИСПРАВЛЕНО: Сценарий 2. Проверка на конкурирующие гипотезы (разница < 1.5)
        elif len(ranked_genes) > 1 and (top_candidate['score'] - ranked_genes.iloc[1]['score']) < 1.5:
            second_candidate = ranked_genes.iloc[1]
            print(f"ВНИМАНИЕ: Неопределенность! Найдено два близких кандидата.")
            print(f"1) {top_candidate['gene']} (Балл: {top_candidate['score']})")
            print(f"   Объяснение: {top_candidate['explanation']}")
            print(f"2) {second_candidate['gene']} (Балл: {second_candidate['score']})")
            
            # Записываем оба варианта
            prediction_str = f"{top_candidate['variant']},{second_candidate['variant']}"
            predictions.append(f"{case_id}\t{prediction_str}\tКонкурирующие гипотезы")
            
        else:
            print(f"УВЕРЕННЫЙ КАНДИДАТ: {top_candidate['gene']} ({top_candidate['variant']}) | Балл: {top_candidate['score']}")
            print(f"Объяснение: {top_candidate['explanation']}") # ДОБАВЛЕНО: Вывод объяснения
            predictions.append(f"{case_id}\t{top_candidate['variant']}\tУверенный кандидат")

    # 5. СОХРАНЕНИЕ РЕЗУЛЬТАТОВ
    output_file = "predictions.tsv"
    with open(output_file, "w") as f:
        # ИСПРАВЛЕНО: Добавлена колонка explanation для жюри
        f.write("case_id\tprediction\treason\n") 
        for line in predictions:
            f.write(f"{line}\n")
    
    print(f"\n{'='*50}")
    print(f"ГОТОВО! Все результаты сохранены в {output_file}")
    print(f"{'='*50}")

if __name__ == "__main__":
    run_pipeline()
