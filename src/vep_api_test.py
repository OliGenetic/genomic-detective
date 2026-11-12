import requests

def get_gene_for_variant(chrom, pos, ref, alt):
    # Адрес API Ensembl (версия GRCh38)
    url = f"https://rest.ensembl.org/vep/human/region/{chrom}:{pos}-{pos}/{alt}?assembly=GRCh38"
    
    headers = {"Content-Type": "application/json"}
    
    print(f"Спрашиваем у Ensembl про вариант: {chrom}:{pos}:{ref}:{alt}...")
    response = requests.get(url, headers=headers)
    
    if response.status_code == 200:
        data = response.json()
        if not data:
            print("VEP ничего не нашел (возможно, это не мутация, а совпадение с референсом).")
            return None
        
        # Достаем название гена из ответа
        result = data[0]
        gene = "Неизвестно"
        if "transcript_consequences" in result:
            for t in result["transcript_consequences"]:
                if "gene_symbol" in t:
                    gene = t["gene_symbol"]
                    break
                    
        print(f"УСПЕХ! Этот вариант находится в гене: {gene}")
        return gene
    else:
        print(f"Ошибка API: {response.status_code}")
        return None

if __name__ == "__main__":
    # Берем вариант, который мы точно знаем из файла с ответами (case_01)
    # chr15:48412568:C:A
    get_gene_for_variant("15", "48412568", "C", "A")