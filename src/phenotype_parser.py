import json

# Путь к файлу с симптомами (JSON)
file_path = "data/case_01_phenotype.json"

def parse_phenotype(path):
    with open(path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    positive_hpo = set()
    negative_hpo = set()
    
    # Проходим по всем симптомам пациента
    for feature in data.get("phenotypicFeatures", []):
        hpo_id = feature["type"]["id"]
        if feature.get("excluded", False):
            negative_hpo.add(hpo_id)
        else:
            positive_hpo.add(hpo_id)
            
    print(f"Найдено положительных симптомов (есть у пациента): {len(positive_hpo)}")
    print(f"Найдено исключенных симптомов (нет у пациента): {len(negative_hpo)}")
    return positive_hpo, negative_hpo

if __name__ == "__main__":
    pos, neg = parse_phenotype(file_path)
    print("\nПримеры симптомов, которые ЕСТЬ:")
    print(list(pos)[:3]) # Покажем первые 3