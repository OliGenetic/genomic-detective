"""
Патч: дозапрашивает варианты, которые остались NA (или NaN) после основного прогона.
Использование:
    python fix_na_variants.py data/case_03_annotated.tsv [data/case_05_annotated.tsv ...]
"""
import pandas as pd
import requests
import time
import sys
import os
from concurrent.futures import ThreadPoolExecutor, as_completed

SESSION = requests.Session()
SESSION.headers.update({"Content-Type": "application/json"})


def annotate_with_retry(chrom, pos, ref, alt, max_retries=5):
    chrom = str(chrom).replace('chr', '')
    pos_int = int(float(pos))
    ref_u = str(ref).strip().upper()
    alt_u = str(alt).strip().upper()

    for a in [alt_u, ref_u]:
        if not a or a in ("NA", ".", "-"):
            continue
        url = f"https://rest.ensembl.org/vep/human/region/{chrom}:{pos_int}-{pos_int}/{a}?assembly=GRCh38&check_ref=0"

        for attempt in range(max_retries):
            try:
                resp = SESSION.get(url, timeout=60)
                if resp.status_code == 200:
                    data = resp.json()
                    if data:
                        item = data[0]
                        gene = "NA"
                        if "transcript_consequences" in item:
                            for t in item["transcript_consequences"]:
                                if "gene_symbol" in t:
                                    gene = t["gene_symbol"]
                                    break
                        if gene and gene != "NA":
                            return gene, item.get("most_severe_consequence", "NA")
                    return "NA", "NA"
                elif resp.status_code == 429:
                    wait = 3 * (2 ** attempt)
                    time.sleep(wait)
                    continue
                elif resp.status_code >= 500:
                    time.sleep(3)
                    continue
                else:
                    break
            except requests.exceptions.Timeout:
                time.sleep(2)
                continue
            except Exception:
                break

    return "NA", "NA"


def patch_file(path):
    print(f"\n{'='*60}")
    print(f">>> ПАТЧ: {path}")
    print(f"{'='*60}")

    if not os.path.exists(path):
        print(f"  ФАЙЛ НЕ НАЙДЕН: {path}")
        return

    df = pd.read_csv(path, sep="\t")
    if 'gene' not in df.columns:
        print("  ОШИБКА: нет колонки 'gene'")
        return

    # ИСПРАВЛЕНО: NaN тоже считаем как NA
    df['gene'] = df['gene'].fillna('NA')
    missing = df[df['gene'] == 'NA'].index.tolist()
    print(f"  NA вариантов: {len(missing)} из {len(df)}")

    if not missing:
        print("  Нечего чинить — всё аннотировано")
        return

    rescued = 0
    rescued_list = []

    def _fetch(idx):
        row = df.loc[idx]
        g, c = annotate_with_retry(row['chrom'], row['pos'], row['ref'], row['alt'])
        return idx, g, c

    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = {pool.submit(_fetch, idx): idx for idx in missing}
        for k, fut in enumerate(as_completed(futures, timeout=2400)):
            try:
                idx, g, c = fut.result(timeout=10)
                if g != "NA":
                    df.at[idx, 'gene'] = g
                    df.at[idx, 'consequence'] = c
                    rescued += 1
                    rescued_list.append((idx, g))
            except Exception:
                pass

            if (k + 1) % 20 == 0:
                print(f"  ... {k+1}/{len(missing)} (спасено: {rescued})")

    df.to_csv(path, sep="\t", index=False)
    print(f"  ✅ Спасено: {rescued}/{len(missing)}")

    if rescued_list:
        genes_found = set(g for _, g in rescued_list)
        print(f"  Уникальных генов среди спасённых: {len(genes_found)}")

        # Спецпроверки для конкретных кейсов
        if "case_03" in path and "KMT2D" in genes_found:
            print(f"  🎯 KMT2D НАЙДЕН! Готов к скорингу.")
        if "case_05" in path and "LMNA" in genes_found:
            print(f"  🎯 LMNA НАЙДЕН! Готов к скорингу.")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Использование: python fix_na_variants.py data/case_XX_annotated.tsv [...]")
        sys.exit(1)

    for path in sys.argv[1:]:
        patch_file(path)
