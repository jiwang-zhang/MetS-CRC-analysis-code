#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Deep verification of GEO downloads in 19.databases.
Checks: file presence, platform id, expression-matrix presence & numeric-ness,
sample counts, and RNA-seq count-file integrity.
"""
import os, gzip, re, sys

BASE = r"D:\daixie-colorectal\19.databases"

# expected: gse -> (expected_platform_gpl, type, expected_samples)
EXPECT = {
    "GSE98895":  ("GPL570",  "array",  40),
    "GSE200744": ("GPL28038","rnaseq", 100),
    "GSE145412": ("GPL18573","rnaseq", 32),
    "GSE21510":  ("GPL570",  "array",  148),
    "GSE32323":  ("GPL570",  "array",  34),
    "GSE39582":  ("GPL570",  "array",  585),
    "GSE10950":  ("GPL6104", "array",  48),
    "GSE21815":  ("GPL6480", "array",  141),
}

def gz_lines(path):
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            yield line.rstrip("\n")

def parse_series_matrix(path):
    """Return dict with platform, sample_count, has_table, n_probes, probe_fmt, n_col_error."""
    platform = None
    sample_ids = []
    has_table = False
    n_probes = 0
    probe_fmt = None
    header_cols = 0
    numeric_ok = True
    table_begin = False
    rows_checked = 0
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("!Platform_ID") or line.startswith("!Sample_platform_id"):
                m = re.search(r'"(GPL\d+)"', line)
                if m: platform = m.group(1)
            if line.startswith("!Series_sample_id"):
                # e.g. !Series_sample_id	"GSM1 GSM2 ..."
                parts = line.split("\t")
                if len(parts) > 1:
                    sample_ids = re.findall(r'GSM\d+', parts[1])
            if line.startswith("!series_matrix_table_begin"):
                table_begin = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if table_begin:
                if n_probes == 0:
                    # header row
                    cols = line.split("\t")
                    header_cols = len(cols) - 1  # first col = ID_REF
                else:
                    cols = line.split("\t")
                    if probe_fmt is None and cols:
                        pid = cols[0]
                        if pid.startswith("ILMN_"): probe_fmt = "ILMN_* (Illumina)"
                        elif re.match(r'^\d+$', pid): probe_fmt = "numeric (Affy/Agilent)"
                        elif pid.startswith(("A_","GE")): probe_fmt = "vendor probe id"
                        else: probe_fmt = pid[:20]
                    # check numeric-ness of a few value cols
                    vals = cols[1:]
                    if rows_checked < 5 and vals:
                        nonnum = [v for v in vals if not re.match(r'^-?\d+(\.\d+)?([eE]-?\d+)?$', v.strip())]
                        if nonnum:
                            numeric_ok = False
                    rows_checked += 1
                n_probes += 1
    return {
        "platform": platform,
        "sample_ids_in_meta": len(sample_ids),
        "has_table": table_begin,
        "n_probes": n_probes,
        "header_sample_cols": header_cols,
        "probe_fmt": probe_fmt,
        "numeric_ok": numeric_ok,
    }

def parse_counts(path, max_rows=2000):
    """Check an RNA-seq count file: shape + integer-valued."""
    n_genes = 0
    n_samples = 0
    all_int = True
    first = True
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for i, line in enumerate(f):
            line = line.rstrip("\n")
            if not line: continue
            cols = line.split("\t")
            if first:
                n_samples = len(cols) - 1
                first = False
                continue
            n_genes += 1
            if n_genes <= 50:  # sample first 50 genes for integer check
                vals = cols[1:]
                for v in vals:
                    v = v.strip()
                    if v == "": continue
                    try:
                        iv = float(v)
                        if abs(iv - round(iv)) > 1e-9:
                            all_int = False
                    except ValueError:
                        all_int = False
            if n_genes >= max_rows:
                break
    return {"n_genes_sampled": min(n_genes, max_rows), "n_samples": n_samples, "all_integer": all_int}

def main():
    print("=" * 78)
    print("GEO DOWNLOAD VERIFICATION  (19.databases)")
    print("=" * 78)
    problems = []
    for gse, (exp_plat, typ, exp_samp) in EXPECT.items():
        d = os.path.join(BASE, gse)
        print(f"\n### {gse}  (expected {typ}, {exp_plat}, ~{exp_samp} samples)")
        sm = os.path.join(d, f"{gse}_series_matrix.txt.gz")
        if not os.path.exists(sm):
            print("  [ERROR] series_matrix.txt.gz MISSING")
            problems.append(f"{gse}: series_matrix missing")
            continue
        info = parse_series_matrix(sm)
        print(f"  series_matrix platform : {info['platform']}  (expected {exp_plat})")
        print(f"  sample_ids in metadata : {info['sample_ids_in_meta']}")
        print(f"  expression table       : {'YES' if info['has_table'] else 'NO (metadata-only)'}")
        print(f"  probe/gene rows        : {info['n_probes']}")
        if info['has_table']:
            print(f"  sample cols in table   : {info['header_sample_cols']}")
            print(f"  probe id format        : {info['probe_fmt']}")
            print(f"  values numeric         : {info['numeric_ok']}")
        # platform mismatch
        if info['platform'] and info['platform'] != exp_plat:
            problems.append(f"{gse}: platform {info['platform']} != expected {exp_plat}")
            print(f"  [WARN] PLATFORM MISMATCH")
        # RNA-seq: check counts file
        if typ == "rnaseq":
            cands = [f for f in os.listdir(d) if 'count' in f.lower() or 'rnaseq' in f.lower()]
            cands = [c for c in cands if c.endswith('.gz') and 'series_matrix' not in c]
            if not cands:
                print("  [ERROR] no RNA-seq count file found")
                problems.append(f"{gse}: missing count file")
                continue
            cf = os.path.join(d, cands[0])
            ci = parse_counts(cf)
            print(f"  count file             : {cands[0]}")
            print(f"  count samples (cols-1): {ci['n_samples']}")
            print(f"  count genes sampled    : {ci['n_genes_sampled']}")
            print(f"  all values integer     : {ci['all_integer']}")
            if not ci['all_integer']:
                problems.append(f"{gse}: count values not integer")
            if ci['n_samples'] != exp_samp:
                problems.append(f"{gse}: count sample {ci['n_samples']} != expected {exp_samp}")
        else:
            # array: verify GPL annot present
            annots = [f for f in os.listdir(d) if 'annot' in f.lower() and f.endswith('.gz')]
            print(f"  GPL annot files        : {annots}")
            if not annots:
                problems.append(f"{gse}: no GPL annot file")
        # sample count check (meta)
        if info['sample_ids_in_meta'] and info['sample_ids_in_meta'] != exp_samp:
            problems.append(f"{gse}: meta sample {info['sample_ids_in_meta']} != expected {exp_samp}")
    print("\n" + "=" * 78)
    if problems:
        print(f"PROBLEMS FOUND ({len(problems)}):")
        for p in problems:
            print("  - " + p)
    else:
        print("NO STRUCTURAL PROBLEMS DETECTED.")
    print("=" * 78)

if __name__ == "__main__":
    main()
