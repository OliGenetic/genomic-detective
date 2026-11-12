import pandas as pd
import requests
import time

def annotate_variants_batch(input_path, output_path, batch_size=200):
    print(f"Загружаем данные из {input_path}...")
    df = pd.read_csv(input_path, sep="\t")
    
    url = "https://rest.ensembl.org/vep/human/region"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    
    genes = []
    consequences = []
    
    for i in range(0, len(df), batch_size):
        batch = df.iloc[i:i+batch_size]
        
        # Формируем список для VEP: "chr:pos:ref:alt"
        # VEP ожидает хромосому без префикса "chr"
        variants_list = []
        for _, row in batch.iterrows():
            chrom = str(row['chrom']).replace('chr', '')
            variants_list.append(f"{chrom}:{row['pos']}:{row['ref']}:{row['alt']}")
        
        payload = {"variants": variants_list}
        print(f"Отправляем батч {i//batch_size + 1} ({len(variants_list)} вариантов)...")
        
        try:
            response = requests.post(url, headers=headers, json=payload, timeout=120)
            
            if response.status_code == 200:
                data = response.json()
                # VEP возвращает результаты в том же порядке
                # Но если вариант не найден, он может его пропустить, поэтому создаем словарь
                results_dict = {}
                for item in data:
                    variant_id = item.get("input", "")
                    # Достаем самый важный тип мутации
                    most_severe = item.get("most_severe_consequence", "NA")
                    
                    gene_symbol = "NA"
                    if "transcript_consequences" in item:
                        for transcript in item["transcript_consequences"]:
                            if "gene_symbol" in transcript:
                                gene_symbol = transcript["gene_symbol"]
                                break
                    
                    # Сохраняем в словарь по позиции
                    results_dict[variant_id] = (gene_symbol, most_severe)
                
                # Сопоставляем с нашим батчем
                for _, row in batch.iterrows():
                    chrom = str(row['chrom']).replace('chr', '')
                    var_id = f"{chrom}:{row['pos']}:{row['ref']}:{row['alt']}"
                    
                    if var_id in results_dict:
                        g, c = results_dict[var_id]
                        genes.append(g)
                        consequences.append(c)
                    else:
                        genes.append("NA")
                        consequences.append("NA")
            else:
                print(f"Ошибка API: {response.status_code} - {response.text}")
                genes.extend(["NA"] * len(batch))
                consequences.extend(["NA"] * len(batch))
                
        except Exception as e:
            print(f"Сетевая ошибка: {e}")
            genes.extend(["NA"] * len(batch))
            consequences.extend(["NA"] * len(batch))
        
        # Пауза, чтобы не забанили API
        time.sleep(2)
    
    # Добавляем колонки и сохраняем
    df['gene'] = genes
    df['consequence'] = consequences
    
    print(f"Сохранение результата в {output_path}...")
    df.to_csv(output_path, sep="\t", index=False)
    print("Готово!")
    
    # Показываем, сколько генов нашли
    found_genes = df[df['gene'] != 'NA']
    print(f"Аннотировано генов: {len(found_genes)} из {len(df)}")
    print("Примеры найденных генов:")
    print(found_genes[['chrom', 'pos', 'gene', 'consequence']].head(10))

if __name__ == "__main__":
    # Входной и выходной файлы
    input_file = "data/case_01_variants.tsv"
    output_file = "data/case_01_annotated.tsv"
    
    annotate_variants_batch(input_file, output_file)