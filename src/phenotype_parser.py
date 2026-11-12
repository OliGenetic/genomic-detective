import json

def parse_phenotype(path):
    with open(path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    positive_hpo = set()
    negative_hpo = set()
    positive_labels = []
    negative_labels = []
    
    for feature in data.get("phenotypicFeatures", []):
        hpo_id = feature["type"]["id"]
        label = feature["type"].get("label", hpo_id) # Берем название, если есть
        
        if feature.get("excluded", False):
            negative_hpo.add(hpo_id)
            negative_labels.append(label)
        else:
            positive_hpo.add(hpo_id)
            positive_labels.append(label)
            
    return positive_hpo, negative_hpo, positive_labels, negative_labels

if __name__ == "__main__":
    pos, neg, pos_l, neg_l = parse_phenotype("data/case_01_phenotype.json")
    print("Есть:", pos_l[:3])
    print("Исключено:", neg_l[:3])