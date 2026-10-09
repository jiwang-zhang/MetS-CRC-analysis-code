# -*- coding: utf-8 -*-
"""对 enrichr ORA 命中的 Top 通路, 计算其重叠 co-DEG 在 CRC / MetS 两病里的
平均 log2FC -> 判断该通路整体"上调"还是"下调"方向。"""
import os, pandas as pd, numpy as np

BASE = r"G:/BaiduSyncdisk/daixie-colorectal/第二部分DEG"
OUT  = os.path.join(BASE, "15.enrich_200coDEG")
csv  = os.path.join(BASE, "09.两个疾病分别DEG", "co_DEGs_CRC_MetS.csv")

df = pd.read_csv(csv).set_index("Gene")
df.index = df.index.str.upper()

def direction(row):
    genes = [g.strip().upper() for g in str(row["Genes"]).split(";") if g.strip()]
    sub = df.reindex([g for g in genes if g in df.index])
    if sub.empty:
        return np.nan, np.nan, ""
    mc = sub["CRC_log2FC"].mean()
    mm = sub["MetS_log2FC"].mean()
    both = sub[(sub["CRC_log2FC"]>0)==(sub["MetS_log2FC"]>0)]  # 同向基因数
    return round(mc,3), round(mm,3), f"{len(both)}/{len(sub)}"

for name, fname in [("GO", "enrichr_GO_top25.csv"), ("KEGG", "enrichr_KEGG_top25.csv")]:
    p = os.path.join(OUT, fname)
    if not os.path.exists(p):
        print(f"[{name}] 缺失 {p}"); continue
    d = pd.read_csv(p)
    res = []
    for _, r in d.iterrows():
        mc, mm, same = direction(r)
        res.append({"Term": r["Term"], "P": r["P-value"], "AdjP": r["Adjusted P-value"],
                    "Overlap": r["Overlap"], "CRC_meanFC": mc, "MetS_meanFC": mm,
                    "同向基因": same})
    res = pd.DataFrame(res)
    res["整体方向"] = np.where(res["CRC_meanFC"]>0, "↑上调", "↓下调")
    out = os.path.join(OUT, f"direction_{name}_top25.csv")
    res.to_csv(out, index=False)
    print(f"\n========== {name} 方向性 (Top, 按 AdjP) ==========")
    print(res[["Term","整体方向","CRC_meanFC","MetS_meanFC","同向基因","AdjP"]].head(15).to_string(index=False))
    print(f"-> 已存 {out}")
