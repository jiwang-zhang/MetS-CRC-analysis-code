# -*- coding: utf-8 -*-
"""把 hotspot_intersections.csv 中命中的热点集按主题合并为"研究方向候选包"，输出去重基因数+基因清单。"""
import csv, os

BASE = r"D:\daixie-colorectal\第二部分DEG\10.获取研究的基因集"
CSV = os.path.join(BASE, "hotspot_intersections.csv")

rows = []
with open(CSV, encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        rows.append(r)

def theme_of(name):
    up = name.upper()
    if any(k in up for k in ["TNFA","INFLAMMATORY","NFKB","NF_KAPPA","CYTOKINE","IL6","IL2","STAT",
                             "INTERFERON","ALLOGRAFT","COMPLEMENT","IMMUNE","REJECTION","TLR",
                             "T CELL","B CELL","CHEST"]):
        return "A_炎症免疫"
    if any(k in up for k in ["FATTY_ACID","CHOLESTEROL","LIPID","PPAR","ADIPOGENESIS","ADIPO","BILE"]):
        return "B_脂质代谢"
    if any(k in up for k in ["E2F","KRAS","G2M","MYC","MITOTIC","SPINDLE","DNA_REPAIR","P53","G1"]):
        return "C_增殖细胞周期"
    if any(k in up for k in ["APOPTOS","FERROPTOS","PYROPTOS","NECROPTOS","CUPROPTOS","DEATH"]):
        return "D_细胞死亡"
    if any(k in up for k in ["INSULIN","WNT","MTOR","PI3K","AMPK","FOXO"]):
        return "E_胰岛素Wnt代谢"
    if any(k in up for k in ["ESTROGEN","ANDROGEN"]):
        return "F_激素"
    return None

themes = {}
for r in rows:
    t = theme_of(r["gene_set"])
    if t is None:
        continue
    themes.setdefault(t, {"sets": [], "genes": set()})
    themes[t]["sets"].append(r["gene_set"])
    themes[t]["genes"].update(r["hit_genes"].split(";"))

print("=== 按主题合并后的研究方向候选包 ===\n")
for t in sorted(themes):
    info = themes[t]
    print(f"[{t}]  去重候选基因 {len(info['genes'])} 个 | 来自 {len(info['sets'])} 个热点集")
    print(f"   基因: {', '.join(sorted(info['genes']))}")
    print()

# 保存每个主题的候选基因文件
for t in sorted(themes):
    info = themes[t]
    out = os.path.join(BASE, f"candidate_{t}.txt")
    with open(out, "w", encoding="utf-8") as f:
        f.write("\n".join(sorted(info["genes"])))
    print(f"  已存 {out}")

# 额外：单集最大命中(TOP1) 与 经典炎症单集(TNFA_NFKB) 单独列出供对比
print("\n=== 单独对比单集 ===")
for target in ["REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM","HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_FATTY_ACID_METABOLISM"]:
    for r in rows:
        if r["gene_set"] == target:
            print(f"{target}: {r['n_hits']} 基因 -> {r['hit_genes']}")
