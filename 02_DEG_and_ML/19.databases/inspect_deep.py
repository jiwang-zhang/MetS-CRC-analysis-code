#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Deep qualitative check on the ambiguous datasets."""
import os, gzip, re

BASE = r"D:\daixie-colorectal\19.databases"

def head_lines(path, n=5):
    out=[]
    with gzip.open(path,"rt",encoding="utf-8",errors="replace") as f:
        for i,line in enumerate(f):
            if i>=n: break
            out.append(line.rstrip("\n"))
    return out

def series_platform_and_samp(path):
    plat=None; samps=0; table_begin=False; nrow=0; ncol=0; pid=None
    with gzip.open(path,"rt",encoding="utf-8",errors="replace") as f:
        for line in f:
            if line.startswith("!Platform_ID") or line.startswith("!Sample_platform_id"):
                m=re.search(r'"(GPL\d+)"',line); 
                if m: plat=m.group(1)
            if line.startswith("!Series_sample_id"):
                samps=len(re.findall(r'GSM\d+',line))
            if line.startswith("!series_matrix_table_begin"): table_begin=True; continue
            if line.startswith("!series_matrix_table_end"): break
            if table_begin:
                cols=line.rstrip("\n").split("\t")
                if nrow==0: ncol=len(cols)-1
                else:
                    if pid is None: pid=cols[0]
                nrow+=1
    return plat,samps,table_begin,nrow,ncol,pid

print("===== GSE98895 series_matrix actual platform =====")
p,s,t,nr,nc,pid=series_platform_and_samp(os.path.join(BASE,"GSE98895","GSE98895_series_matrix.txt.gz"))
print(f"  declared platform = {p}")
print(f"  sample_ids in meta = {s}")
print(f"  has expr table     = {t}, probes={nr}, sample_cols={nc}, probe_id_example={pid}")
print("  note: GSE98895 is Illumina HumanHT-12 -> GPL6947. If 'p' above is GPL6947, the")
print("        download was CORRECT; the leftover GPL570.annot.gz is just an unused wrong file.")

print("\n===== GSE200744 count file: real value distribution =====")
with gzip.open(os.path.join(BASE,"GSE200744","GSE200744_RNAseq_count_in_Mets.txt.gz"),"rt") as f:
    header=f.readline().rstrip("\n").split("\t")
    ncol=len(header)
    # sample first 3 data rows fully
    rows=[]
    for _ in range(3):
        rows.append(f.readline().rstrip("\n").split("\t"))
vals=[]
for r in rows:
    for v in r[1:]:
        try: vals.append(float(v))
        except: pass
import statistics
nonint=sum(1 for v in vals if abs(v-round(v))>1e-9)
print(f"  columns = {ncol} (1 id + {ncol-1} samples)")
print(f"  header[0..3] = {header[:4]}")
print(f"  sampled {len(vals)} values; min={min(vals):.4g} max={max(vals):.4g} mean={statistics.mean(vals):.4g}")
print(f"  non-integer values among sample = {nonint}/{len(vals)}")
print("  -> if nonint>0 and values are fractional (e.g. 0.xx..), this is a NORMALIZED matrix (TPM/FPKM), NOT raw counts.")

print("\n===== GSE145412 count file: real value distribution =====")
with gzip.open(os.path.join(BASE,"GSE145412","GSE145412_rnaseq.txt.gz"),"rt") as f:
    header=f.readline().rstrip("\n").split("\t")
    ncol=len(header)
    rows=[]
    for _ in range(3):
        rows.append(f.readline().rstrip("\n").split("\t"))
vals=[]
for r in rows:
    for v in r[1:]:
        try: vals.append(float(v))
        except: pass
nonint=sum(1 for v in vals if abs(v-round(v))>1e-9)
print(f"  columns = {ncol} (1 id + {ncol-1} samples)")
print(f"  header[0..3] = {header[:4]}")
print(f"  sampled {len(vals)} values; min={min(vals):.4g} max={max(vals):.4g} mean={statistics.mean(vals):.4g}")
print(f"  non-integer values among sample = {nonint}/{len(vals)}")
