import pandas as pd
import os

def calculate_gene_scores(annotated_tsv_path, positive_hpo, negative_hpo, hpo_db_path="data/genes_to_phenotype.txt"):
    if not os.path.exists(hpo_db_path):
        print(f"ОШИБКА: Файл базы HPO не найден по пути {hpo_db_path}")
        return None

    print("Загружаем базу HPO... (это может занять пару секунд)")
    hpo_db = pd.read_csv(hpo_db_path, sep="\t", comment="#")
    gene_to_hpo = hpo_db.groupby('gene_symbol')['hpo_id'].apply(set).to_dict()
    
    print("Загружаем аннотированные варианты...")
    df = pd.read_csv(annotated_tsv_path, sep="\t")
    df_genes = df[df['gene'] != 'NA'].copy()
    
    results = []
    for _, row in df_genes.iterrows():
        gene = row['gene']
        consequence = str(row['consequence'])
        
        # --- 1. БИОЛОГИЧЕСКИЙ ФИЛЬТР ПО ТИПУ МУТАЦИИ ---
        cons_bonus = 0
        if consequence in ['frameshift_variant', 'stop_gained', 'splice_acceptor_variant', 'splice_donor_variant']:
            cons_bonus = 3.0  # Обрыв белка - очень вероятно, что это причина
        elif consequence in ['missense_variant', 'inframe_deletion', 'inframe_insertion']:
            cons_bonus = 2.0  # Замена аминокислоты - частый виновник
        elif consequence in ['synonymous_variant', 'intron_variant', 'intergenic_variant', 'upstream_gene_variant', 'downstream_gene_variant']:
            cons_bonus = -2.0 # Тихие мутации - скорее всего, это просто шум

        # --- 1.5 БОНУС ЗА CLINVAR (Если вариант уже известен врачам) ---
        # Проверяем, есть ли колонка клинической значимости в таблице
        clinvar_sig = str(row.get('clinvar_significance', '')).lower()
        if 'pathogenic' in clinvar_sig and 'conflicting' not in clinvar_sig:
            cons_bonus += 5.0 # Огромный бонус, если это точно патогенный вариант

        # --- 2. ФЕНОТИПИЧЕСКИЙ ФИЛЬТР (СИМПТОМЫ) ---
        if gene not in gene_to_hpo:
            score = 0
            matches = 0
            excluded_matches = 0
        else:
            gene_hpo_set = gene_to_hpo[gene]
            matches = len(gene_hpo_set.intersection(positive_hpo))
            excluded_matches = len(gene_hpo_set.intersection(negative_hpo))
            
            # Формула: Совпадения*1.5 + Бонус_за_мутацию - Штраф_за_исключенные
            score = (matches * 1.5) + cons_bonus - (excluded_matches * 0.5)
        
        results.append({
            'gene': gene,
            'variant': f"{row['chrom']}:{row['pos']}:{row['ref']}:{row['alt']}",
            'consequence': consequence,
            'score': round(score, 2),
            'matches': matches,
            'excluded_matches': excluded_matches,
            'cons_bonus': cons_bonus
        })
    
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values(by='score', ascending=False)
    return results_df

if __name__ == "__main__":
    import sys
    sys.path.append('src')
    from phenotype_parser import parse_phenotype
    
    print("--- ЗАПУСК ОБНОВЛЕННОГО АЛГОРИТМА СКОРИНГА (v2.1 с ClinVar) ---")
    pos_hpo, neg_hpo = parse_phenotype("data/case_01_phenotype.json")
    
    ranked_genes = calculate_gene_scores(
        annotated_tsv_path="data/case_01_annotated.tsv",
        positive_hpo=pos_hpo,
        negative_hpo=neg_hpo
    )
    
    if ranked_genes is not None and not ranked_genes.empty:
        print("\n--- ТОП-5 ГЕНОВ-КАНДИДАТОВ (case_01) ---")
        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', 1000)
        print(ranked_genes.head(5).to_string(index=False))