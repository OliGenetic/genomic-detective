import pandas as pd
import requests
import time
import os

def get_clinvar_sig_ncbi(chrom, pos, ref, alt):
    """
    Запрашивает клиническую значимость в NCBI ClinVar через E-utilities.
    Возвращает строку ('Pathogenic', 'Benign', 'Unknown').
    """
    clean_chrom = str(chrom).replace('chr', '')
    search_term = f"{clean_chrom}:{pos} {ref}>{alt}"
    
    # Шаг 1: Ищем ID варианта в ClinVar
    search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": "clinvar",
        "term": search_term,
        "retmode": "json"
    }
    
    try:
        response = requests.get(search_url, params=params, timeout=10)
        if response.status_code != 200:
            return "Unknown"
            
        data = response.json()
        id_list = data.get("esearchresult", {}).get("idlist", [])
        
        if not id_list:
            return "Unknown" # NCBI не знает такой вариант
            
        variation_id = id_list[0]
        
        # Шаг 2: Запрашиваем сводку по ID
        summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        sum_params = {
            "db": "clinvar",
            "id": variation_id,
            "retmode": "json"
        }
        
        sum_response = requests.get(summary_url, params=sum_params, timeout=10)
        if sum_response.status_code == 200:
            sum_data = sum_response.json()
            result = sum_data.get("result", {}).get(variation_id, {})
            germline = result.get("germline_classification", {})
            return germline.get("description", "Unknown")
            
    except Exception:
        pass
        
    return "Unknown"


def annotate_variants_batch(input_path, output_path, batch_size=200):
    print(f"Загружаем данные из {input_path}...")
    if not os.path.exists(input_path):
        print(f"ОШИБКА: Файл {input_path} не найден.")
        return

    df = pd.read_csv(input_path, sep="\t")
    
    url = "https://rest.ensembl.org/vep/human/region"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    
    genes = []
    consequences = []
    
    # --- ЭТАП 1: АННОТАЦИЯ ЧЕРЕЗ ENSEMBL VEP ---
    for i in range(0, len(df), batch_size):
        batch = df.iloc[i:i+batch_size]
        
        variants_list = []
        for _, row in batch.iterrows():
            chrom = str(row['chrom']).replace('chr', '')
            variants_list.append(f"{chrom}:{row['pos']}:{row['ref']}:{row['alt']}")
        
        payload = {"variants": variants_list}
        print(f"Отправляем батч {i//batch_size + 1} ({len(variants_list)} вариантов) в VEP...")
        
        try:
            response = requests.post(url, headers=headers, json=payload, timeout=120)
            
            if response.status_code == 200:
                data = response.json()
                results_dict = {}
                for item in data:
                    variant_id = item.get("input", "")
                    most_severe = item.get("most_severe_consequence", "NA")
                    
                    gene_symbol = "NA"
                    if "transcript_consequences" in item:
                        for transcript in item["transcript_consequences"]:
                            if "gene_symbol" in transcript:
                                gene_symbol = transcript["gene_symbol"]
                                break
                    
                    results_dict[variant_id] = (gene_symbol, most_severe)
                
                for _, row in batch.iterrows():
                    chrom = str(row['chrom']).replace('chr', '')
                    var_id = f"{chrom}:{row['pos']}:{row['ref']}:{row['alt']}"
                    
                    found = False
                    for key, (g, c) in results_dict.items():
                        if str(row['pos']) in key and row['ref'] in key and row['alt'] in key:
                            genes.append(g)
                            consequences.append(c)
                            found = True
                            break
                    
                    if not found:
                        genes.append("NA")
                        consequences.append("NA")
            else:
                print(f"Ошибка VEP API: {response.status_code}")
                genes.extend(["NA"] * len(batch))
                consequences.extend(["NA"] * len(batch))
                
        except Exception as e:
            print(f"Сетевая ошибка VEP: {e}")
            genes.extend(["NA"] * len(batch))
            consequences.extend(["NA"] * len(batch))
        
        time.sleep(3) # Пауза для VEP
    
    df['gene'] = genes
    df['consequence'] = consequences

    # --- ЭТАП 2: ПРОВЕРКА CLINVAR ЧЕРЕЗ NCBI ---
    print("\nНачинаем проверку клинической значимости в NCBI ClinVar (это займет время)...")
    clinvar_significances = []
    
    for index, row in df.iterrows():
        if row['gene'] != 'NA':
            sig = get_clinvar_sig_ncbi(row['chrom'], row['pos'], row['ref'], row['alt'])
            clinvar_significances.append(sig)
            # Задержка 0.4 сек, чтобы NCBI не выдал ошибку 429 (Too Many Requests)
            time.sleep(0.4)
        else:
            clinvar_significances.append("NA")
            
    df['clinvar_significance'] = clinvar_significances

    print(f"Сохранение результата в {output_path}...")
    df.to_csv(output_path, sep="\t", index=False)
    print("Готово!")
    
    found_genes = df[df['gene'] != 'NA']
    print(f"Аннотировано генов: {len(found_genes)} из {len(df)}")
    
    path_vars = df[df['clinvar_significance'].str.contains('athogenic', case=False, na=False)]
    if not path_vars.empty:
         print(f"Найдено патогенных вариантов по ClinVar: {len(path_vars)}")
    else:
         print("Патогенных вариантов по ClinVar не найдено (или API не ответил).")


if __name__ == "__main__":
    input_file = "data/case_01_variants.tsv"
    output_file = "data/case_01_annotated.tsv"
    
    annotate_variants_batch(input_file, output_file)