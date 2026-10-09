# -*- coding: utf-8 -*-
"""
将 MSigDB v2025.1 全库(xlsm) 中 MetS×CRC 相关热点基因集 与 200 个 co-DEG 取交集。
输出：每个热点集的命中基因数、命中基因、命中率，按命中数降序。
原则 (对齐 IJBM 8.5 模板)：用权威下载集(MSigDB)，不自制；热点方向=炎症/胰岛素-IGF/Wnt/脂质代谢/细胞死亡/免疫。
"""
import openpyxl
import csv, os

BASE = r"D:\daixie-colorectal\第二部分DEG"
CO_DEG_CSV = os.path.join(BASE, r"09.两个疾病分别DEG", "co_DEGs_CRC_MetS.csv")
XLSX = os.path.join(BASE, r"10.获取研究的基因集", "如果做其它基因集，可以看这个文件.xlsm")
OUT_CSV = os.path.join(BASE, r"10.获取研究的基因集", "hotspot_intersections.csv")

# ---------- 1. 读 200 co-DEG ----------
co_genes = []
with open(CO_DEG_CSV, encoding="utf-8") as f:
    r = csv.reader(f)
    header = next(r)
    for row in r:
        if row and row[0].strip():
            co_genes.append(row[0].strip().upper())
co_set = set(co_genes)
print(f"[co-DEG] 共 {len(co_set)} 个基因")

# ---------- 2. 定义要保留的热点集(按 collection + 关键字) ----------
# 我们只对与 MetS×CRC 机制相关的集做交集，避免 3.5 万集全部扫描(也全部交了没意义)
KEEP_HALLMARK = True  # Hallmark 50 个全要(它们本身就是"热点签名")

# KEGG / REACTOME 精选关键字(大小写不敏感匹配集名)
KEGG_KEEP = [
    "COLORECTAL CANCER", "INSULIN SIGNALING", "PPAR SIGNALING", "ADIPOCYTOKINE",
    "WNT SIGNALING", "NF KAPPA", "FERROPTOSIS", "APOPTOSIS", "PANCREATIC CANCER",
    "GLYCOLYSIS", "PI3K", "MTOR", "FOXO", "AMPK", "FATTY ACID", "CHOLESTEROL",
    "BILE ACID", "INSULIN RESISTANCE", "TYPE II DIABETES", "ADIPOSE", "INFLAMMATORY",
    "TNF SIGNALING", "IL17", "TH17", "PD L1", "ANTIGEN PROCESSING", "T CELL",
    "B CELL", "NOD LIKE", "TOLL LIKE", "CYTOKINE", "JAK STAT",
]
REACTOME_KEEP = [
    "INSULIN", "TNF", "WNT", "PPAR", "ADIPOGENESIS", "ADIPOCYTE", "INFLAMMATORY",
    "CYTOKINE", "INTERFERON", "T CELL", "B CELL", "INNATE IMMUNE", "ADAPTIVE IMMUNE",
    "METABOLISM OF LIPIDS", "FATTY ACID", "CHOLESTEROL", "SIGNALING BY", "APOPTOSIS",
    "CELL DEATH", "NECROPTOSIS", "PYROPTOSIS", "FERROPTOSIS", "NLRP3", "NF-KB",
    "MTOR", "PI3K", "KRAS", "COLORECTAL", "CANCER", "METABOLISM OF GLUCOSE",
]
GO_KEEP = [
    "INSULIN", "ADIPOGENESIS", "ADIPOCYTE", "FERROPTOSIS", "PYROPTOSIS", "NECROPTOSIS",
    "CUPROPTOSIS", "WNT", "INFLAMMATORY RESPONSE", "TNF", "APOPTOTIC", "LIPID METABOLIC",
    "COLORECTAL", "INTESTINAL", "GLUCOSE", "INSULIN RESISTANCE", "TYPE 2 DIABETES",
    "T CELL", "B CELL", "MACROPHAGE", "CYTOKINE", "INTERFERON", "NF-KAPPAB",
]

def collect_ok(name):
    up = name.upper()
    if up.startswith("HALLMARK_"):
        return KEEP_HALLMARK
    if up.startswith("KEGG_"):
        return any(k in up for k in KEGG_KEEP)
    if up.startswith("REACTOME_"):
        return any(k in up for k in REACTOME_KEEP)
    if up.startswith("GO_"):
        return any(k in up for k in GO_KEEP)
    return False

# ---------- 3. 解析 xlsm(只读流)，取集 ----------
print("[parse] 读取 xlsm(35134 集)...")
wb = openpyxl.load_workbook(XLSX, read_only=True)
ws = wb["msigdb.v2025.1.Hs.symbols"]

sets = {}          # name -> set(genes)
collection_count = {}
for i, row in enumerate(ws.iter_rows(values_only=True), 1):
    if not row or row[0] is None:
        continue
    name = str(row[0]).strip()
    if not collect_ok(name):
        continue
    # 基因从 col3 起 (1-based col index 3 -> row[2])
    genes = {str(g).strip().upper() for g in row[2:] if g is not None and str(g).strip()}
    sets[name] = genes
    coll = name.split("_")[0]
    collection_count[coll] = collection_count.get(coll, 0) + 1
    if i % 5000 == 0:
        print(f"  scanned {i} rows, kept {len(sets)} sets")
wb.close()
print(f"[parse] 保留热点集 {len(sets)} 个: {collection_count}")

# ---------- 4. 交集 ----------
results = []
for name, genes in sets.items():
    hit = sorted(co_set & genes)
    if hit:
        results.append({
            "gene_set": name,
            "collection": name.split("_")[0],
            "n_set_genes": len(genes),
            "n_hits": len(hit),
            "hit_fraction": round(len(hit)/len(genes), 4),
            "hit_genes": ";".join(hit),
        })
results.sort(key=lambda x: (-x["n_hits"], -x["hit_fraction"]))

# ---------- 5. 输出 ----------
with open(OUT_CSV, "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=["gene_set","collection","n_set_genes","n_hits","hit_fraction","hit_genes"])
    w.writeheader()
    for r in results:
        w.writerow(r)

print(f"\n[done] 命中≥1 的热点集共 {len(results)} 个，已存 {OUT_CSV}\n")
print("=== TOP 40 (按命中数) ===")
print(f"{'gene_set':52s} {'coll':10s} {'setN':>5s} {'hit':>4s} {'frac':>6s}  genes")
for r in results[:40]:
    genes_short = r["hit_genes"][:80]
    print(f"{r['gene_set'][:50]:52s} {r['collection']:10s} {r['n_set_genes']:5d} {r['n_hits']:4d} {r['hit_fraction']:6.2f}  {genes_short}")
