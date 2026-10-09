"""
GSEA: hub-gene high/low expression group KEGG enrichment
Python equivalent of OE21.GSEA.R (clusterProfiler GSEA)
Uses gseapy.prerank for GSEA on ranked gene list
"""
import os
import pandas as pd
import numpy as np
import gseapy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================================================
# Paths
# ============================================================
BASE = r"D:\BaiduSyncdisk\daixie-colorectal"
EXP_FILE  = os.path.join(BASE, "第二部分DEG", "25.diagnostic_model", "merged_after_batch_removal.csv")
GENE_FILE = os.path.join(BASE, "4.免疫浸润单细胞", "21.GSEA", "hubGenes.List.txt")
GMT_FILE  = os.path.join(BASE, "4.免疫浸润单细胞", "21.GSEA", "c2.cp.kegg.Hs.symbols.gmt")
OUT_DIR   = os.path.join(BASE, "4.免疫浸润单细胞", "21.GSEA")

os.chdir(OUT_DIR)
print(f"Working dir: {OUT_DIR}")

# ============================================================
# 1. Read expression matrix
# ============================================================
print("Reading expression matrix...")
rt = pd.read_csv(EXP_FILE, index_col=0, header=0)
print(f"  Raw: {rt.shape[0]} genes x {rt.shape[1]} samples")

# Filter tumor samples (_tre suffix)
tre_cols = [c for c in rt.columns if c.endswith("_tre")]
exp = rt[tre_cols].copy()
print(f"  Tumor only: {exp.shape[0]} genes x {exp.shape[1]} samples")

# avereps: average duplicate gene symbols
exp = exp.groupby(exp.index).mean()
print(f"  After avereps: {exp.shape[0]} genes")

# Remove genes with rowMeans == 0
exp = exp[exp.mean(axis=1) > 0]
print(f"  After filtering zero-mean: {exp.shape[0]} genes")

# ============================================================
# 2. Read hub genes
# ============================================================
hub_genes = [line.strip() for line in open(GENE_FILE) if line.strip()]
print(f"Hub genes: {', '.join(hub_genes)}")

# ============================================================
# 3. GSEA for each hub gene
# ============================================================
for gene in hub_genes:
    if gene not in exp.index:
        print(f"Skip {gene}: not in expression matrix")
        continue

    print(f"\n--- Processing {gene} ---")

    # Split samples by median expression
    gene_exp = exp.loc[gene]
    median_val = gene_exp.median()
    low_mask  = gene_exp < median_val
    high_mask = gene_exp >= median_val

    low_data  = exp.loc[:, low_mask]
    high_data = exp.loc[:, high_mask]
    print(f"  Low group: {low_data.shape[1]} samples, High group: {high_data.shape[1]} samples")

    # Compute mean per gene in each group
    mean_low  = low_data.mean(axis=1)
    mean_high = high_data.mean(axis=1)

    # Avoid division by zero
    mean_low = mean_low.clip(lower=0.00001)
    mean_high = mean_high.clip(lower=0.00001)

    # logFC = meanH - meanL (same as R script: difference of means on log2 scale)
    logfc = mean_high - mean_low
    logfc = logfc.sort_values(ascending=False)
    print(f"  Ranked list: {len(logfc)} genes, range [{logfc.min():.4f}, {logfc.max():.4f}]")

    # Save ranked list
    rnk_df = pd.DataFrame({"gene": logfc.index, "score": logfc.values})
    rnk_path = os.path.join(OUT_DIR, f"{gene}.rnk")
    rnk_df.to_csv(rnk_path, sep="\t", index=False, header=False)

    # Run GSEA
    print(f"  Running GSEA with KEGG gene sets...")
    gsea_res = gseapy.prerank(
        rnk=logfc,
        gene_sets=GMT_FILE,
        min_size=5,
        max_size=500,
        permutation_num=1000,
        seed=42,
        threads=4,
        no_plot=True,
        verbose=False,
    )

    # Get results
    res_df = gsea_res.res2d.copy()
    print(f"  GSEA raw results: {len(res_df)} pathways")

    # Filter by pvalue < 0.05 (matching R script)
    # gseapy columns: Term, ES, NES, NOM p-val, FDR q-val, ...
    # R clusterProfiler uses 'pvalue' for NOM p-val
    pval_col = None
    for col in res_df.columns:
        if "p-val" in col.lower() and "fdr" not in col.lower() and "fwer" not in col.lower():
            pval_col = col
            break
    if pval_col is None:
        pval_col = "NOM p-val" if "NOM p-val" in res_df.columns else res_df.columns[3]

    nes_col = "NES" if "NES" in res_df.columns else "NES"

    print(f"  Using pvalue column: '{pval_col}', NES column: '{nes_col}'")

    sig_df = res_df[res_df[pval_col] < 0.05].copy()
    print(f"  Significant (p<0.05): {len(sig_df)} pathways")

    # Save full result
    full_out = os.path.join(OUT_DIR, f"{gene}.result.txt")
    res_df.to_csv(full_out, sep="\t", index=False)
    print(f"  Saved: {gene}.result.txt")

    # Save significant result
    if len(sig_df) > 0:
        sig_out = os.path.join(OUT_DIR, f"{gene}.sig.txt")
        sig_df.to_csv(sig_out, sep="\t", index=False)

    # Plot top 5 upregulated (NES > 0) and top 5 downregulated (NES < 0)
    for direction, label, nes_filter in [
        ("high", f"Enriched in high {gene} group", lambda x: x > 0),
        ("low",  f"Enriched in low {gene} group",  lambda x: x < 0),
    ]:
        subset = sig_df[sig_df[nes_col].apply(nes_filter)] if len(sig_df) > 0 else pd.DataFrame()
        if len(subset) == 0:
            # If no significant, use all results
            subset = res_df[res_df[nes_col].apply(nes_filter)].head(5)
        if len(subset) == 0:
            print(f"  No {direction} pathways to plot")
            continue

        # Sort by absolute NES
        subset = subset.reindex(subset[nes_col].abs().sort_values(ascending=False).index)
        top_terms = subset.head(5)

        pdf_path = os.path.join(OUT_DIR, f"{gene}.{direction}Exp.pdf")

        fig, axes = plt.subplots(len(top_terms), 1, figsize=(6.5, 2.5 * len(top_terms)),
                                 constrained_layout=True)
        if len(top_terms) == 1:
            axes = [axes]

        for idx, (_, row) in enumerate(top_terms.iterrows()):
            term = row["Term"]
            term_key = term

            # Get enrichment plot data from gseapy
            try:
                # gseapy stores results in .results as dict
                if term in gsea_res.results:
                    res = gsea_res.results[term]
                    hit_idx = res["hits"]
                    es_vals = res["RES"]
                    ranked = logfc.values

                    x = np.arange(len(es_vals))
                    axes[idx].plot(x, es_vals, color="#2c7fb8", linewidth=1.2)
                    axes[idx].fill_between(x, es_vals, 0, alpha=0.2, color="#2c7fb8")

                    # Mark hit positions
                    for h in hit_idx:
                        if h < len(x):
                            axes[idx].axvline(x=h, color="#d95f0e", alpha=0.3, linewidth=0.5)

                    nes_val = row[nes_col]
                    pval_val = row[pval_col]
                    short_term = term if len(term) < 60 else term[:57] + "..."
                    axes[idx].set_title(f"{short_term}\nNES={nes_val:.2f}, p={pval_val:.3g}",
                                        fontsize=7)
                    axes[idx].set_ylabel("Enrichment Score", fontsize=7)
                    axes[idx].tick_params(labelsize=6)
                else:
                    axes[idx].text(0.5, 0.5, f"{term}\n(plot data not available)",
                                  ha="center", va="center", fontsize=8)
                    axes[idx].set_xticks([])
                    axes[idx].set_yticks([])
            except Exception as e:
                axes[idx].text(0.5, 0.5, f"{term}\nError: {e}",
                              ha="center", va="center", fontsize=7)
                axes[idx].set_xticks([])
                axes[idx].set_yticks([])

        fig.suptitle(label, fontsize=9, fontweight="bold")
        fig.savefig(pdf_path, format="pdf", dpi=300)
        plt.close(fig)
        print(f"  Saved: {gene}.{direction}Exp.pdf ({len(top_terms)} terms)")

    print(f"  {gene} GSEA done")

print("\n=== GSEA all done ===")
