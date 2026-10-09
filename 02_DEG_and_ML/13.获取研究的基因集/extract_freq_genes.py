# ============================================================
# 从 HALLMARK_TNFA_SIGNALING_VIA_NFKB 创始集 GMT 中
# 提取出现 ≥3 次的高频基因
# ============================================================
from collections import Counter

import os
os.chdir(r"D:\daixie-colorectal\第二部分DEG\10.获取研究的基因集")
gmt_file = "HALLMARK_TNFA_SIGNALING_VIA_NFKB_FOUNDERS.v2026.1.Hs.gmt"
min_freq  = 3  # 可改：2/3/5/10

# 读 GMT，统计每个基因在几个创始集中出现
gene_freq = Counter()
total_sets = 0
with open(gmt_file) as f:
    for line in f:
        parts = line.strip().split('\t')
        name  = parts[0]
        genes = set(parts[2:])  # 同一创始集内去重
        for g in genes:
            gene_freq[g] += 1
        total_sets += 1

print(f"创始集总数: {total_sets}")
print(f"去重基因总数: {len(gene_freq)}")

# 筛选 ≥ min_freq 的基因
selected = {g for g, c in gene_freq.items() if c >= min_freq}
print(f"出现 ≥{min_freq} 次的基因: {len(selected)}")

# 保存为 TXT（每行一个基因）
out_txt = f"TNF_NFKB_freq{min_freq}.txt"
with open(out_txt, "w") as f:
    for g in sorted(selected):
        f.write(g + "\n")
print(f"已保存: {out_txt}")

# 同时保存带频次的详细版
out_csv = f"TNF_NFKB_freq{min_freq}_with_counts.csv"
with open(out_csv, "w") as f:
    f.write("Gene,FounderSetCount\n")
    for g, c in sorted(gene_freq.items(), key=lambda x: -x[1]):
        if c >= min_freq:
            f.write(f"{g},{c}\n")
print(f"已保存: {out_csv}")

# 统计频次分布
freq_dist = Counter(gene_freq.values())
print(f"\n频次分布:")
for fq in sorted(freq_dist, reverse=True)[:10]:
    print(f"  出现 {fq} 次: {freq_dist[fq]} 个基因")
