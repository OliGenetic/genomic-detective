import pandas as pd
import requests
import time
import os
import sys
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed

# =========================================================================
# КОНСТАНТЫ
# =========================================================================
VEP_POST_URL = "https://rest.ensembl.org/vep/human/region"
VEP_GET_URL = "https://rest.ensembl.org/vep/human/region/{chrom}:{pos}-{pos}/{alt}"
CLINVAR_SEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
CLINVAR_SUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}

SESSION = requests.Session()
SESSION.headers.update(HEADERS)

CLINVAR_CACHE = {}


# =========================================================================
# CLINVAR
# =========================================================================
def get_clinvar_sig(chrom, pos, ref, alt):
    key = f"{chrom}:{pos}:{ref}:{alt}"
    if key in CLINVAR_CACHE:
        return CLINVAR_CACHE[key]

    clean_chrom = str(chrom).replace('chr', '')
    search_term = f"{clean_chrom}:{pos} {ref}>{alt}"
    params = {"db": "clinvar", "term": search_term, "retmode": "json"}

    try:
        r = SESSION.get(CLINVAR_SEARCH, params=params, timeout=15)
        if r.status_code != 200:
            CLINVAR_CACHE[key] = "Unknown"
            return "Unknown"

        id_list = r.json().get("esearchresult", {}).get("idlist", [])
        if not id_list:
            CLINVAR_CACHE[key] = "Unknown"
            return "Unknown"

        vid = id_list[0]
        r2 = SESSION.get(
            CLINVAR_SUMMARY,
            params={"db": "clinvar", "id": vid, "retmode": "json"},
            timeout=15,
        )
        if r2.status_code == 200:
            result = r2.json().get("result", {}).get(vid, {})
            sig = result.get("germline_classification", {}).get("description", "Unknown")
            CLINVAR_CACHE[key] = sig
            return sig
    except Exception:
        pass

    CLINVAR_CACHE[key] = "Unknown"
    return "Unknown"


# =========================================================================
# VEP: GET (5 retry, timeout=60, экспоненциальный backoff)
# =========================================================================
def annotate_single_get(chrom, pos, ref, alt, max_retries=5):
    chrom = str(chrom).replace('chr', '')
    pos_int = int(float(pos))
    ref_u = str(ref).strip().upper()
    alt_u = str(alt).strip().upper()

    for a in [alt_u, ref_u]:
        if not a or a in ("NA", ".", "-"):
            continue
        url = VEP_GET_URL.format(chrom=chrom, pos=pos_int, alt=a) + "?assembly=GRCh38&check_ref=0"

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
                    print(f"      [429] ждём {wait} сек (попытка {attempt+1}/{max_retries})")
                    time.sleep(wait)
                    continue
                elif resp.status_code >= 500:
                    time.sleep(3)
                    continue
                else:
                    break
            except requests.exceptions.Timeout:
                print(f"      [timeout] {url[-40:]} (попытка {attempt+1}/{max_retries})")
                time.sleep(2)
                continue
            except Exception:
                break

    return "NA", "NA"


# =========================================================================
# ПАРСИНГ POST-БАТЧА
# =========================================================================
def _parse_vep_response(data):
    out = {}
    for item in data:
        try:
            input_str = str(item.get("input", "")).strip()
            if input_str and ":" in input_str:
                parts = input_str.split(":")
                if len(parts) == 4:
                    pos = parts[1].strip()
                    ref = parts[2].strip().upper()
                    alt = parts[3].strip().upper()

                    gene = "NA"
                    if "transcript_consequences" in item:
                        for t in item["transcript_consequences"]:
                            if "gene_symbol" in t:
                                gene = t["gene_symbol"]
                                break

                    cons = item.get("most_severe_consequence", "NA")
                    out[(pos, ref, alt)] = (gene, cons)
                    continue

            pos = str(int(item.get("start")))
            allele_str = str(item.get("allele_string", "")).upper()
            if "/" in allele_str:
                ref = allele_str.split("/")[0]
                alts = allele_str.split("/")[1].split(",")

                gene = "NA"
                if "transcript_consequences" in item:
                    for t in item["transcript_consequences"]:
                        if "gene_symbol" in t:
                            gene = t["gene_symbol"]
                            break

                cons = item.get("most_severe_consequence", "NA")
                for alt in alts:
                    out[(pos, ref, alt.strip())] = (gene, cons)
        except Exception:
            continue
    return out


# =========================================================================
# ОСНОВНАЯ ФУНКЦИЯ
# =========================================================================
def annotate_variants_batch(input_path, output_path, batch_size=100, use_clinvar=True):
    t0 = time.time()
    print(f"\n{'='*60}")
    print(f">>> АННОТАЦИЯ: {input_path}")
    print(f"{'='*60}")

    if not os.path.exists(input_path):
        print(f"ОШИБКА: Файл {input_path} не найден.")
        return

    df = pd.read_csv(input_path, sep="\t")
    if 'genotype' not in df.columns:
        df['genotype'] = "Unknown"

    n = len(df)
    print(f"    Вариантов в файле: {n}")

    genes = ["NA"] * n
    consequences = ["NA"] * n
    missing_indices = []

    # --- ЭТАП 1: POST-батчи ---
    n_batches = (n + batch_size - 1) // batch_size
    for bi in range(n_batches):
        start = bi * batch_size
        end = min(start + batch_size, n)

        variants_list = []
        for idx in range(start, end):
            row = df.iloc[idx]
            chrom = str(row['chrom']).replace('chr', '')
            pos = str(int(float(row['pos'])))
            ref = str(row['ref']).strip().upper()
            alt = str(row['alt']).strip().upper()
            variants_list.append(f"{chrom}:{pos}:{ref}:{alt}")

        try:
            resp = SESSION.post(VEP_POST_URL, json={"variants": variants_list}, timeout=120)
            if resp.status_code != 200:
                print(f"    [Батч {bi+1}/{n_batches}] ОШИБКА {resp.status_code}: {resp.text[:150]}")
                missing_indices.extend(range(start, end))
                time.sleep(1)
                continue

            parsed = _parse_vep_response(resp.json())

            batch_found = 0
            for idx in range(start, end):
                row = df.iloc[idx]
                pos = str(int(float(row['pos'])))
                ref = str(row['ref']).strip().upper()
                alt = str(row['alt']).strip().upper()

                if (pos, ref, alt) in parsed:
                    genes[idx], consequences[idx] = parsed[(pos, ref, alt)]
                    batch_found += 1
                elif (pos, alt, ref) in parsed:
                    genes[idx], consequences[idx] = parsed[(pos, alt, ref)]
                    batch_found += 1
                else:
                    missing_indices.append(idx)

            print(f"    [Батч {bi+1}/{n_batches}] OK | Найдено: {batch_found}/{end-start}")
        except Exception as e:
            print(f"    [Батч {bi+1}/{n_batches}] Exception: {e}")
            missing_indices.extend(range(start, end))

        time.sleep(1)

    # --- ЭТАП 2: ДВА РАУНДА GET ---
    if missing_indices:
        def _fetch_one(idx):
            row = df.iloc[idx]
            g, c = annotate_single_get(row['chrom'], row['pos'], row['ref'], row['alt'])
            return idx, g, c

        def _run_round(indices, max_workers, label):
            rescued_local = 0
            with ThreadPoolExecutor(max_workers=max_workers) as pool:
                futures = {pool.submit(_fetch_one, idx): idx for idx in indices}
                try:
                    for k, fut in enumerate(as_completed(futures, timeout=2400)):
                        try:
                            idx, g, c = fut.result(timeout=5)
                            if g != "NA" and genes[idx] == "NA":
                                genes[idx] = g
                                consequences[idx] = c
                                rescued_local += 1
                        except Exception:
                            pass
                        if (k + 1) % 50 == 0:
                            print(f"      [{label}] {k+1}/{len(indices)} (спасено: {rescued_local})")
                except TimeoutError:
                    print(f"      [{label}] общий таймаут")
            return rescued_local

        print(f"\n    Раунд 1: GET для {len(missing_indices)} вариантов (6 потоков)...")
        rescued = _run_round(missing_indices, max_workers=6, label="GET-1")
        print(f"    Раунд 1: спасено {rescued}/{len(missing_indices)}")

        still_missing = [i for i in missing_indices if genes[i] == "NA"]
        if still_missing:
            print(f"\n    Раунд 2: GET для {len(still_missing)} упорных (3 потока, timeout=60)...")
            rescued2 = _run_round(still_missing, max_workers=3, label="GET-2")
            print(f"    Раунд 2: спасено {rescued2}/{len(still_missing)}")
            rescued += rescued2

        print(f"    Итого спасено через GET: {rescued}/{len(missing_indices)}")

    # --- ЭТАП 3: ClinVar ---
    if use_clinvar:
        to_check = [i for i in range(n) if genes[i] != "NA"]
        print(f"\n    Проверка ClinVar для {len(to_check)} генов (6 потоков)...")

        def _fetch_clinvar(idx):
            row = df.iloc[idx]
            sig = get_clinvar_sig(row['chrom'], row['pos'], row['ref'], row['alt'])
            return idx, sig

        clinvar_sig = ["NA"] * n
        with ThreadPoolExecutor(max_workers=6) as pool:
            futures = {pool.submit(_fetch_clinvar, idx): idx for idx in to_check}
            try:
                for k, fut in enumerate(as_completed(futures, timeout=2400)):
                    try:
                        idx, sig = fut.result(timeout=5)
                        clinvar_sig[idx] = sig
                    except Exception:
                        pass
                    if (k + 1) % 50 == 0:
                        print(f"      ... проверено {k+1}/{len(to_check)} (кэш: {len(CLINVAR_CACHE)})")
            except TimeoutError:
                print("      !!! Таймаут ClinVar")

        df['clinvar_significance'] = clinvar_sig
    else:
        df['clinvar_significance'] = "Unknown"
        print("\n    [ClinVar пропущен — быстрый режим]")

    df['gene'] = genes
    df['consequence'] = consequences

    # --- ИСПРАВЛЕНО: защита от NaN перед сохранением ---
    df['gene'] = df['gene'].where(df['gene'].notna() & (df['gene'] != ''), 'NA')
    df['consequence'] = df['consequence'].where(df['consequence'].notna() & (df['consequence'] != ''), 'NA')

    df.to_csv(output_path, sep="\t", index=False)

    found = sum(1 for g in genes if g != "NA")
    path_vars = sum(
        1 for s in df['clinvar_significance']
        if 'athogenic' in str(s).lower() and 'conflicting' not in str(s).lower()
    )
    elapsed = time.time() - t0

    print(f"\n    ✅ ИТОГО за {elapsed:.1f} сек:")
    print(f"       Аннотировано генов: {found}/{n}")
    print(f"       Патогенных по ClinVar: {path_vars}")
    print(f"       Сохранено: {output_path}")


# =========================================================================
# CLI
# =========================================================================
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--no-clinvar", action="store_true",
                        help="Пропустить ClinVar (быстрее в 3-5 раз)")
    args = parser.parse_args()
    annotate_variants_batch(args.input, args.output, use_clinvar=not args.no_clinvar)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        main()
    else:
        annotate_variants_batch("data/case_02_variants.tsv", "data/case_02_annotated.tsv")
