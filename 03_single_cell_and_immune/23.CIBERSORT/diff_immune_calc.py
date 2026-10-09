import csv, statistics, sys
from scipy import stats
from statsmodels.stats.multitest import multipletests

path = r"D:\BaiduSyncdisk\daixie-colorectal\4.免疫浸润单细胞\23.CIBERSORT\CIBERSORT-Results.txt"

with open(path, newline="") as f:
    reader = csv.reader(f, delimiter="\t")
    header = next(reader)
    rows = [r for r in reader if r and r[0].strip()]

cells = header[1:]  # 22 immune cells
print("n_cols(header)=", len(header), "n_cells=", len(cells), "n_rows=", len(rows))

# group by suffix
normal, tumor = [], []
for r in rows:
    sid = r[0].strip()
    vals = [float(x) for x in r[1:]]
    (normal if sid.endswith("_con") else tumor).append(vals)

print(f"Normal(_con)={len(normal)}  Tumor(_tre)={len(tumor)}  total={len(normal)+len(tumor)}")

results = []
for i, cell in enumerate(cells):
    nv = [row[i] for row in normal]
    tv = [row[i] for row in tumor]
    try:
        U, p = stats.mannwhitneyu(tv, nv, alternative="two-sided")
    except ValueError:
        p = float("nan")
    med_n, med_t = statistics.median(nv), statistics.median(tv)
    mean_n, mean_t = statistics.mean(nv), statistics.mean(tv)
    direction = "CRC↑" if mean_t > mean_n else "CRC↓"
    results.append((cell, med_n, med_t, mean_n, mean_t, mean_t - mean_n, p))

pvals = [r[6] for r in results]
rej, qvals, _, _ = multipletests(pvals, method="fdr_bh")

print(f"\nSignificant at p<0.05: {sum(1 for p in pvals if p<0.05)} / 22")
print(f"Significant at BH FDR<0.05: {sum(rej)} / 22\n")

print(f"{'Cell':32s} {'mean_Nor':>8s} {'mean_CRC':>8s} {'dir':>5s} {'p':>10s} {'q(BH)':>10s}  sig")
print("-"*95)
for (cell, mn, mt, mean_n, mean_t, d, p), q, rj in sorted(zip(results, qvals, rej), key=lambda x: x[0][6]):
    flag = "***" if p<1e-3 else ("**" if p<1e-2 else ("*" if p<0.05 else ""))
    sig = "BH+" if rj else ""
    dirs = "CRC↑" if mean_t > mean_n else "CRC↓"
    print(f"{cell:32s} {mean_n:8.4f} {mean_t:8.4f} {dirs:>5s} {p:10.2e} {q:10.2e}  {flag}{sig}")

# focused: significant cells grouped by direction (mean-based)
sig_up = [r[0] for r, q, rj in zip(results, qvals, rej) if rj and r[5] > 0]
sig_down = [r[0] for r, q, rj in zip(results, qvals, rej) if rj and r[5] < 0]
print("\n--- BH-significant, CRC-enriched (mean) ---")
print(", ".join(sig_up))
print("--- BH-significant, CRC-depleted (mean) ---")
print(", ".join(sig_down))
