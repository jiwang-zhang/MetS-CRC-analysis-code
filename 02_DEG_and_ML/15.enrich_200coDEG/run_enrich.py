# -*- coding: utf-8 -*-
"""
用 200 个 co-DEGs (MetS ∩ CRC 共同差异基因) 跑 GO / KEGG 方向分析。
- 输入: 09 文件夹的 co_DEGs_CRC_MetS.csv (Gene, CRC_log2FC, CRC_adjP, MetS_log2FC, MetS_adjP)
- 方法1: gseapy.enrichr  ORA (直接等价于 clusterProfiler enrichGO/enrichKEGG；需联网 Enrichr)
- 方法2: gseapy.gsea + 本地 gmt (离线可靠；带方向, 用 signed 排序)
排序指标: score = sign(mean_log2FC) * ( -log10(CRC_adjP) - log10(MetS_adjP) ) / 2
  -> 正值=两病共同显著且偏向上调; 负值=偏向下调。GSEA 据此分上下尾给方向。
"""
import os, sys
import pandas as pd
import numpy as np

BASE = r"G:/BaiduSyncdisk/daixie-colorectal/第二部分DEG"
CSV  = os.path.join(BASE, "09.两个疾病分别DEG", "co_DEGs_CRC_MetS.csv")
GO_GMT  = os.path.join(BASE, "12.GSEA", "c5.go.Hs.symbols.gmt")
KEGG_GMT = os.path.join(BASE, "12.GSEA", "c2.cp.kegg.Hs.symbols.gmt")
OUT  = os.path.join(BASE, "15.enrich_200coDEG")
os.makedirs(OUT, exist_ok=True)

df = pd.read_csv(CSV)
print(f"[info] 读入 {len(df)} 个 co-DEGs")
genes = df["Gene"].astype(str).tolist()

# ---- 构建带方向的 pre-ranked 列表 (GSEA 用) ----
p_crc  = -np.log10(df["CRC_adjP"].clip(lower=1e-300))
p_mets = -np.log10(df["MetS_adjP"].clip(lower=1e-300))
mean_fc = (df["CRC_log2FC"] + df["MetS_log2FC"]) / 2.0
score = np.sign(mean_fc) * (p_crc + p_mets) / 2.0
rank = pd.DataFrame({"gene": df["Gene"].astype(str).str.upper(), "score": score})
rank = rank.sort_values("score", ascending=False)
rank.to_csv(os.path.join(OUT, "preranked_200coDEGs.txt"), sep="\t", index=False, header=False)
print(f"[info] pre-ranked 列表已写: {len(rank)} 基因 (正={int((score>0).sum())}, 负={int((score<0).sum())})")

import gseapy as gp

# ================= 方法1: enrichr ORA (GO BP / KEGG) =================
print("\n========== [1] enrichr ORA (联网) ==========")
try:
    res_go = gp.enrichr(gene_list=genes,
                         gene_sets=["GO_Biological_Process_2021"],
                         organism="human",
                         outdir=os.path.join(OUT, "enrichr_GO"),
                         no_plot=True, verbose=False)
    res_kegg = gp.enrichr(gene_list=genes,
                           gene_sets=["KEGG_2021_Human"],
                           organism="human",
                           outdir=os.path.join(OUT, "enrichr_KEGG"),
                           no_plot=True, verbose=False)
    go_df = res_go.results
    kegg_df = res_kegg.results
    go_df = go_df.sort_values("Adjusted P-value").head(25)
    kegg_df = kegg_df.sort_values("Adjusted P-value").head(25)
    go_df.to_csv(os.path.join(OUT, "enrichr_GO_top25.csv"), index=False)
    kegg_df.to_csv(os.path.join(OUT, "enrichr_KEGG_top25.csv"), index=False)
    print(f"[OK] enrichr GO: {len(res_go.results)} 条显著; KEGG: {len(res_kegg.results)} 条显著")
    print("\n--- GO Biological Process Top15 (按 Adjusted P-value) ---")
    print(go_df[["Term","Overlap","P-value","Adjusted P-value","Genes"]].head(15).to_string(index=False))
    print("\n--- KEGG Top15 (按 Adjusted P-value) ---")
    print(kegg_df[["Term","Overlap","P-value","Adjusted P-value","Genes"]].head(15).to_string(index=False))
except Exception as e:
    print(f"[WARN] enrichr 联网失败: {e}\n        -> 改用本地 GSEA (方法2)")

# ================= 方法2: 本地 gmt GSEA (带方向, 离线) =================
print("\n========== [2] 本地 gmt GSEA (离线, 带方向) ==========")
for name, gmt in [("GO", GO_GMT), ("KEGG", KEGG_GMT)]:
    try:
        pre = os.path.join(OUT, "preranked_200coDEGs.txt")
        gs = gp.prerank(rnk=pre, gene_sets=gmt,
                         permutation_num=1000, method="permutation",
                         outdir=os.path.join(OUT, f"gsea_{name}"), no_plot=True, verbose=False)
        res = gs.res2d
        res = res.sort_values("fdr", ascending=False)  # 先按fdr升序更直观, 下面再排
        res = res.sort_values("fdr")
        # 分开上下尾
        up = res[res["nes"] > 0].sort_values("fdr").head(20)
        dn = res[res["nes"] < 0].sort_values("fdr").head(20)
        up.to_csv(os.path.join(OUT, f"gsea_{name}_UP_top20.csv"), index=False)
        dn.to_csv(os.path.join(OUT, f"gsea_{name}_DOWN_top20.csv"), index=False)
        print(f"\n--- {name} GSEA: 共 {len(res)} 通路, 显著(fdr<0.25)={int((res['fdr']<0.25).sum())} ---")
        print(f"\n  >> {name} 上调方向 Top12 (NES>0):")
        print(up[["term","nes","pval","fdr","matched_size","ledge_genes"]].head(12).to_string(index=False))
        print(f"\n  >> {name} 下调方向 Top12 (NES<0):")
        print(dn[["term","nes","pval","fdr","matched_size","ledge_genes"]].head(12).to_string(index=False))
    except Exception as e:
        print(f"[ERR] {name} GSEA 失败: {e}")

print("\n[done] 所有结果已写入:", OUT)
