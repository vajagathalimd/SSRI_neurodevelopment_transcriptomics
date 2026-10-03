#!/usr/bin/env python3

import csv
import json
import time
import urllib.error
import urllib.request
from pathlib import Path

PROJECT = Path("/mnt/d/SSRI_neurodevelopment_transcriptomics")

DE_DIR = PROJECT / "09_differential_expression" / "GSE166297"
OUT_DIR = PROJECT / "11_pathway_enrichment" / "GSE166297_gprofiler"
OUT_DIR.mkdir(parents=True, exist_ok=True)

API_URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"

def read_deg_csv(path):
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)
    return rows

def clean_gene_id(g):
    if g is None:
        return ""
    g = str(g).strip()
    if g.upper() in ["NA", "NAN", "NULL", ""]:
        return ""
    if g.startswith("ENSG") and "." in g:
        g = g.split(".")[0]
    return g

def get_gene(row):
    gene_name = clean_gene_id(row.get("gene_name", ""))
    gene_id = clean_gene_id(row.get("gene_id", ""))

    if gene_name:
        return gene_name
    return gene_id

def unique_genes(rows):
    genes = []
    seen = set()

    for row in rows:
        g = get_gene(row)

        if not g:
            continue

        # Avoid problematic non-gene labels
        if g.lower().startswith("novel"):
            continue

        if g not in seen:
            genes.append(g)
            seen.add(g)

    return genes

def split_by_direction(rows):
    up = []
    down = []

    for row in rows:
        try:
            lfc = float(row.get("log2FoldChange", "nan"))
        except Exception:
            continue

        if lfc > 0:
            up.append(row)
        elif lfc < 0:
            down.append(row)

    return up, down

def write_gene_list(genes, name):
    out_file = OUT_DIR / f"{name}_input_genes.txt"
    with open(out_file, "w", encoding="utf-8") as f:
        for g in genes:
            f.write(g + "\n")

def call_gprofiler(genes, name):
    print(f"Running g:Profiler for {name}: {len(genes)} genes")
    write_gene_list(genes, name)

    if len(genes) < 5:
        print(f"Skipping {name}: fewer than 5 genes")
        return []

    payload = {
        "organism": "hsapiens",
        "query": genes,
        "sources": ["GO:BP", "GO:MF", "GO:CC", "REAC", "KEGG", "WP"],
        "user_threshold": 0.05,
        "significance_threshold_method": "g_SCS",
        "all_results": False,
        "ordered": False,
        "no_evidences": False,
        "domain_scope": "annotated"
    }

    data = json.dumps(payload).encode("utf-8")

    req = urllib.request.Request(
        API_URL,
        data=data,
        headers={
            "Content-Type": "application/json",
            "User-Agent": "GSE166297-transcriptomics-analysis"
        },
        method="POST"
    )

    result = None

    for attempt in range(1, 4):
        try:
            with urllib.request.urlopen(req, timeout=180) as response:
                result = json.loads(response.read().decode("utf-8"))
            break

        except urllib.error.HTTPError as e:
            body = e.read().decode("utf-8", errors="replace")
            print(f"Attempt {attempt} failed for {name}: HTTP {e.code}")
            print("Server response:")
            print(body[:2000])
            time.sleep(10)

        except Exception as e:
            print(f"Attempt {attempt} failed for {name}: {e}")
            time.sleep(10)

    if result is None:
        print(f"Failed enrichment for {name}")
        return []

    terms = result.get("result", [])

    if not terms:
        print(f"No significant terms for {name}")
        return []

    out_file = OUT_DIR / f"{name}_gprofiler_enrichment.csv"

    fieldnames = [
        "source",
        "native",
        "term_name",
        "p_value",
        "significant",
        "term_size",
        "query_size",
        "intersection_size",
        "effective_domain_size",
        "precision",
        "recall",
        "parents",
        "intersections"
    ]

    with open(out_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for term in terms:
            intersections = term.get("intersections", [])
            parents = term.get("parents", [])

            writer.writerow({
                "source": term.get("source", ""),
                "native": term.get("native", ""),
                "term_name": term.get("name", ""),
                "p_value": term.get("p_value", ""),
                "significant": term.get("significant", ""),
                "term_size": term.get("term_size", ""),
                "query_size": term.get("query_size", ""),
                "intersection_size": term.get("intersection_size", ""),
                "effective_domain_size": term.get("effective_domain_size", ""),
                "precision": term.get("precision", ""),
                "recall": term.get("recall", ""),
                "parents": json.dumps(parents, ensure_ascii=False) if isinstance(parents, (list, dict)) else parents,
                "intersections": json.dumps(intersections, ensure_ascii=False) if isinstance(intersections, (list, dict)) else intersections
            })

    print(f"Saved: {out_file} | significant terms: {len(terms)}")
    return terms

ic5_file = DE_DIR / "FLX_IC5_vs_Control_D10_significant_FDR05_log2FC058.csv"
ic20_file = DE_DIR / "FLX_IC20_100_vs_Control_D10_significant_FDR05_log2FC058.csv"

ic5_rows = read_deg_csv(ic5_file)
ic20_rows = read_deg_csv(ic20_file)

ic5_genes = unique_genes(ic5_rows)
ic20_genes = unique_genes(ic20_rows)
shared_genes = sorted(set(ic5_genes).intersection(set(ic20_genes)))

ic5_up_rows, ic5_down_rows = split_by_direction(ic5_rows)
ic20_up_rows, ic20_down_rows = split_by_direction(ic20_rows)

gene_sets = {
    "FLX_IC5_all_significant_DEGs": ic5_genes,
    "FLX_IC20_100_all_significant_DEGs": ic20_genes,
    "Shared_FLX_IC5_and_IC20_100_DEGs": shared_genes,
    "FLX_IC5_upregulated_DEGs": unique_genes(ic5_up_rows),
    "FLX_IC5_downregulated_DEGs": unique_genes(ic5_down_rows),
    "FLX_IC20_100_upregulated_DEGs": unique_genes(ic20_up_rows),
    "FLX_IC20_100_downregulated_DEGs": unique_genes(ic20_down_rows)
}

summary_file = OUT_DIR / "GSE166297_gprofiler_input_gene_counts.csv"

with open(summary_file, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    writer.writerow(["gene_set", "n_genes"])
    for name, genes in gene_sets.items():
        writer.writerow([name, len(genes)])

print(f"Saved input summary: {summary_file}")

for name, genes in gene_sets.items():
    call_gprofiler(genes, name)

print("g:Profiler API enrichment completed.")
