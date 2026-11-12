import pandas as pd
import requests
import time
import os

def get_clinvar_sig_ncbi(chrom, pos, ref, alt):
    """Запрос в NCBI ClinVar (как было раньше)"""
    clean_chrom = str(chrom).replace('chr', '')
    search_term = f"{clean_chrom}:{pos} {ref}>{alt}"
    search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {"db": "clinvar", "term": search_term, "retmode": "json"}
    try:
        response = requests.get(search_url, params=params, timeout=10)
        if response.status_code != 200: return "Unknown"
        data = response.json()
        id_list = data.get("esearchresult", {}).get("idlist", [])
        if not id_list: return "Unknown"
        
        variation_id = id_list[0]
        summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        sum_params = {"db": "clinvar", "id": variation_id, "retmode": "json"}
        sum_response = requests.get(summary_url, params=sum_params, timeout=10)
        if sum_response.status_code == 200:
            sum_data = sum_response.json()
            result = sum_data.get("result", {}).get(variation_id, {})
            germline = result.get("germline_classification", {})
            return germline.get("description", "Unknown")
    except Exception:
        pass
    return "Unknown"

def annotate_single_variant(chrom, pos, ref, alt):
    """Индивидуальный запрос для 'потерянных' вариантов"""
    chrom = str(chrom).replace('chr', '')
    url = f"https://rest.ensembl.org/vep/human/region/{chrom}:{pos}-{pos}/{alt}?assembly=GRCh38"
    headers = {"Content-Type": "application/json"}
    try:
        response = requests.get(url, headers=headers, timeout=10)
        if response.status_code == 200:
            data = response.json()
            if data:
                item = data[0]
                most_severe = item.get("most_severe_consequence", "NA")
                gene_symbol = "NA"
                if "transcript_consequences" in item:
                    for transcript in item["transcript_consequences"]:
                        if "gene_symbol" in transcript:
                            gene_symbol = transcript["gene_symbol"]
                            break
                return gene_symbol, most_severe
    except Exception:
        pass
    return "NA", "NA"

def annotate_variants_batch(input_path, output_path, batch_size=100):
    print(f"Загружаем данные из {input_path}...")
    if not os.path.exists(input_path):
        print(f"ОШИБКА: Файл {input_path} не найден.")
        return

    df = pd.read_csv(input_path, sep="\t")
    
    url = "https://rest.ensembl.org/vep/human/region"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    
    genes = []
    consequences = []
    
    # --- ЭТАП 1: АННОТАЦИЯ БАТЧАМИ ---
    for i in range(0, len(df), batch_size):
        batch = df.iloc[i:i+batch_size]
        variants_list = []
        for _, row in batch.iterrows():
            chrom = str(row['chrom']).replace('chr', '')
            variants_list.append(f"{chrom}:{row['pos']}:{row['ref']}:{row['alt']}")
        
        payload = {"variants": variants_list}
        print(f"Отправляем батч {i//batch_size + 1} ({len(variants_list)} вариантов) в VEP...")
        
        batch_genes = {}
        try:
            response = requests.post(url, headers=headers, json=payload, timeout=120)
            if response.status_code == 200:
                data = response.json()
                for item in data:
                    variant_id = item.get("input", "")
                    most_severe = item.get("most_severe_consequence", "NA")
                    gene_symbol = "NA"
                    if "transcript_consequences" in item:
                        for transcript in item["transcript_consequences"]:
                            if "gene_symbol" in transcript:
                                gene_symbol = transcript["gene_symbol"]
                                break
                    # Сохраняем по ключу-позиции
                    batch_genes[str(item.get("start"))] = (gene_symbol, most_severe)
        except Exception as e:
            print(f"Ошибка батча: {e}")
        
        # Сверяем: если варианта нет в ответе — запрашиваем отдельно
        for _, row in batch.iterrows():
            pos_str = str(row['pos'])
            if pos_str in batch_genes:
                g, c = batch_genes[pos_str]
                genes.append(g)
                consequences.append(c)
            else:
                # СПАСАТЕЛЬНЫЙ КРУГ: индивидуальный запрос
                print(f"  -> Вариант {pos_str} потерялся в батче. Запрашиваем отдельно...")
                g, c = annotate_single_variant(row['chrom'], row['pos'], row['ref'], row['alt'])
                genes.append(g)
                consequences.append(c)
                time.sleep(0.2) # Пауза между индивидуальными запросами
                
        time.sleep(5) # Увеличили паузу после батча
    
    df['gene'] = genes
    df['consequence'] = consequences

    # --- ЭТАП 2: CLINVAR ---
    print("\nНачинаем проверку ClinVar (NCBI)...")
    clinvar_significances = []
    for index, row in df.iterrows():
        if row['gene'] != 'NA':
            sig = get_clinvar_sig_ncbi(row['chrom'], row['pos'], row['ref'], row['alt'])
            clinvar_significances.append(sig)
            time.sleep(0.4)
        else:
            clinvar_significances.append("NA")
    df['clinvar_significance'] = clinvar_significances

    df.to_csv(output_path, sep="\t", index=False)
    print("Готово!")
    found_genes = df[df['gene'] != 'NA']
    print(f"Аннотировано генов: {len(found_genes)} из {len(df)}")
    path_vars = df[df['clinvar_significance'].str.contains('athogenic', case=False, na=False)]
    print(f"Найдено патогенных вариантов по ClinVar: {len(path_vars)}")

if __name__ == "__main__":
    annotate_variants_batch("data/case_02_variants.tsv", "data/case_02_annotated.tsv")