import pandas as pd
import os


def generate_explanation(row, pos_hpo, neg_hpo):
    reasons = []

    if row['matches'] > 0:
        reasons.append(f"Совпадает с {row['matches']} симптомами пациента")
    if row['excluded_matches'] > 0:
        reasons.append(f"Есть {row['excluded_matches']} исключенных симптомов (штраф)")

    if 'pathogenic' in str(row['clinvar_sig']).lower():
        reasons.append("Вариант подтвержден как патогенный в ClinVar")

    if row['consequence'] in ['frameshift_variant', 'stop_gained', 'splice_acceptor_variant', 'splice_donor_variant']:
        reasons.append(f"Мутация с потерей функции ({row['consequence']})")
    elif row['consequence'] in ['missense_variant']:
        reasons.append("Миссенс-мутация (замена аминокислоты)")

    if row['genotype'] == '1/1':
        reasons.append("Гомозиготное состояние (вероятно рецессивное заболевание)")
    elif row['genotype'] == '0/1' and row['cons_bonus'] >= 3.0:
        reasons.append("Гетерозиготное состояние")

    if not reasons:
        return "Нет явных доказательств в пользу этого гена."

    return "; ".join(reasons) + "."


def calculate_gene_scores(annotated_tsv_path, positive_hpo, negative_hpo, hpo_db_path="data/genes_to_phenotype.txt"):
    if not os.path.exists(hpo_db_path):
        print(f"ОШИБКА: Файл базы HPO не найден по пути {hpo_db_path}")
        return None

    print("Загружаем базу HPO... (это может занять пару секунд)")
    hpo_db = pd.read_csv(hpo_db_path, sep="\t", comment="#")
    gene_to_hpo = hpo_db.groupby('gene_symbol')['hpo_id'].apply(set).to_dict()

    print("Загружаем аннотированные варианты...")
    df = pd.read_csv(annotated_tsv_path, sep="\t")

    if 'gene' in df.columns:
        df['gene'] = df['gene'].fillna('NA')
    if 'consequence' in df.columns:
        df['consequence'] = df['consequence'].fillna('NA')
    if 'clinvar_significance' in df.columns:
        df['clinvar_significance'] = df['clinvar_significance'].fillna('Unknown')

    df_genes = df[df['gene'] != 'NA'].copy()

    results = []
    for _, row in df_genes.iterrows():
        gene = row['gene']
        consequence = str(row['consequence'])
        genotype = str(row.get('genotype', 'Unknown'))

        # 1. Бонус за тип мутации
        cons_bonus = 0
        if consequence in ['frameshift_variant', 'stop_gained',
                           'splice_acceptor_variant', 'splice_donor_variant',
                           'splice_donor_region_variant', 'splice_acceptor_region_variant',
                           'splice_polypyrimidine_tract_variant']:
            cons_bonus = 3.0
        elif consequence in ['missense_variant', 'inframe_deletion', 'inframe_insertion']:
            cons_bonus = 2.0
        elif consequence in ['synonymous_variant', 'intron_variant', 'intergenic_variant',
                             'upstream_gene_variant', 'downstream_gene_variant']:
            cons_bonus = -2.0

        # 2. Анализ наследования (зиготности)
        if genotype in ['1/1', '1|1']:
            cons_bonus += 0.5
        elif genotype in ['0/1', '1/0', '0|1', '1|0']:
            cons_bonus += 1.0

        # 3. Бонус ClinVar
        clinvar_sig = str(row.get('clinvar_significance', '')).lower()
        if 'pathogenic' in clinvar_sig and 'conflicting' not in clinvar_sig:
            cons_bonus += 5.0

        # 4. Фенотипический фильтр
        if gene not in gene_to_hpo:
            score = 0
            matches = 0
            excluded_matches = 0
        else:
            gene_hpo_set = gene_to_hpo[gene]
            matches = len(gene_hpo_set.intersection(positive_hpo))
            excluded_matches = len(gene_hpo_set.intersection(negative_hpo))
            score = (matches * 1.5) + cons_bonus - (excluded_matches * 0.5)

        results.append({
            'gene': gene,
            'variant': f"{row['chrom']}:{row['pos']}:{row['ref']}:{row['alt']}",
            'consequence': consequence,
            'score': round(score, 2),
            'matches': matches,
            'excluded_matches': excluded_matches,
            'cons_bonus': cons_bonus,
            'genotype': genotype,
            'clinvar_sig': clinvar_sig
        })

    results_df = pd.DataFrame(results)
    if results_df.empty:
        return results_df

    results_df = results_df.sort_values(by='score', ascending=False)

    results_df['explanation'] = results_df.apply(
        lambda r: generate_explanation(r, positive_hpo, negative_hpo), axis=1
    )

    return results_df


if __name__ == "__main__":
    import sys
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    from phenotype_parser import parse_phenotype

    print("--- ЗАПУСК СКОРИНГА ---")
    pos_hpo, neg_hpo, _, _ = parse_phenotype("data/case_01_phenotype.json")

    ranked_genes = calculate_gene_scores(
        annotated_tsv_path="data/case_01_annotated.tsv",
        positive_hpo=pos_hpo,
        negative_hpo=neg_hpo
    )

    if ranked_genes is not None and not ranked_genes.empty:
        print("\n--- ТОП-5 ГЕНОВ-КАНДИДАТОВ ---")
        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', 2000)
        pd.set_option('display.max_colwidth', None)
        print(ranked_genes.head(5).to_string(index=False))
    else:
        print("Кандидаты не найдены.")
