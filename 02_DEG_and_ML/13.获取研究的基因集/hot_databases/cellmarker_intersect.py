# -*- coding: utf-8 -*-
"""CellMarker 2.0 (Cell_marker_All.xlsx) 与 200 co-DEG 取交集：
把我们的共差异基因映射到细胞类型(支持后续 scRNA 定位)；并重点看 37 个炎症免疫候选基因。"""
import openpyxl, os, csv
from collections import defaultdict

BASE = r"D:\daixie-colorectal\第二部分DEG"
XLSX = os.path.join(BASE, r"10.获取研究的基因集\hot_databases", "Cell_marker_All.xlsx")
CO_DEG = os.path.join(BASE, r"09.两个疾病分别DEG", "co_DEGs_CRC_MetS.csv")
INFLAM = os.path.join(BASE, r"10.获取研究的基因集", "candidate_A_炎症免疫.txt")
OUT = os.path.join(BASE, r"10.获取研究的基因集\hot_databases")

# 读 co-DEG
co = []
with open(CO_DEG, encoding="utf-8") as f:
    next(f)
    for line in f:
        g = line.split(",")[0].strip().strip('"')
        if g: co.append(g.upper())
co_set = set(co)

# 读 37 炎症免疫候选
infl = []
with open(INFLAM, encoding="utf-8") as f:
    for line in f:
        g = line.strip().upper()
        if g: infl.append(g)
infl_set = set(infl)

# 解析 CellMarker：gene -> {cell_name set, cell_type set}
gene_cellname = defaultdict(set)
gene_celltype = defaultdict(set)
gene_cancer = defaultdict(set)
n_rows = 0
wb = openpyxl.load_workbook(XLSX, read_only=True)
ws = wb["All"]
hdr = None
for row in ws.iter_rows(values_only=True):
    if hdr is None:
        hdr = row
        idx = {name: i for i, name in enumerate(hdr)}
        i_sym = idx["Symbol"]; i_cn = idx["cell_name"]; i_ct = idx["cell_type"]; i_sp = idx["species"]; i_ca = idx["cancer_type"]
        continue
    n_rows += 1
    sym = row[i_sym]
    if not sym:
        continue
    sym = str(sym).strip().upper()
    gene_cellname[sym].add(str(row[i_cn]))
    gene_celltype[sym].add(str(row[i_ct]))
    gene_cancer[sym].add(str(row[i_ca]))
wb.close()
print(f"[CellMarker] 解析 {n_rows} 行, {len(gene_cellname)} 个基因有标记")

# 交集
def report(geneset, name):
    hit = sorted(geneset & set(gene_cellname.keys()))
    print(f"\n=== {name}: {len(hit)}/{len(geneset)} 个基因是已知 CellMarker 细胞标记 ===")
    # 细胞类型计数
    ct_count = defaultdict(int)
    cn_count = defaultdict(int)
    for g in hit:
        for ct in gene_celltype[g]:
            ct_count[ct] += 1
        for cn in gene_cellname[g]:
            cn_count[cn] += 1
    print(f"  涉及细胞类型(cell_type) Top15: {dict(sorted(ct_count.items(), key=lambda x:-x[1])[:15])}")
    print(f"  涉及具体细胞(cell_name) Top15: {dict(sorted(cn_count.items(), key=lambda x:-x[1])[:15])}")
    # 保存明细
    with open(os.path.join(OUT, f"cellmarker_{name}.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["gene","n_celltypes","cell_types","cancer_types"])
        for g in hit:
            w.writerow([g, len(gene_cellname[g]), ";".join(sorted(gene_cellname[g])), ";".join(sorted(gene_cancer[g]))])
    return hit

report(co_set, "coDEG200")
report(infl_set, "inflammation37")

# 全部 co-DEG 中标记免疫相关细胞的基因(额外洞察)
immune_kw = ["T cell","B cell","Macrophage","NK cell","Dendritic","Neutrophil","Monocyte","Mast","Plasma","Myeloid","Lymphoid","Leukocyte","Granulocyte","Microglia"]
print("\n=== co-DEG 中标记免疫细胞的基因 ===")
for g in sorted(co_set & set(gene_cellname.keys())):
    cns = gene_cellname[g]
    imm = [c for c in cns if any(k.lower() in c.lower() for k in immune_kw)]
    if imm:
        print(f"  {g}: {', '.join(imm[:6])}")
