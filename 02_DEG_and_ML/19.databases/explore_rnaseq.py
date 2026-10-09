import gzip, os, csv, re

BASE = r"D:\daixie-colorectal\19.databases"

def explore_count(gse, fname, n_show=3):
    print("="*70)
    print(f"[{gse}] 计数矩阵: {fname}")
    p = os.path.join(BASE, gse, fname)
    with gzip.open(p, 'rt', errors='ignore') as f:
        r = csv.reader(f, delimiter='\t')
        header = next(r)
        print(f"  列数={len(header)}  前6列={header[:6]}  末3列={header[-3:]}")
        rows = []
        for i, row in enumerate(r):
            rows.append(row)
            if i >= n_show-1:
                break
    ng = sum(1 for _ in gzip.open(p,'rt',errors='ignore'))
    print(f"  总行数(基因数,含header)={ng}")
    print(f"  基因ID示例={[r[0] for r in rows]}")
    print(f"  首行前6值={rows[0][:6]}")
    return header

def parse_series_samples(gse):
    print("="*70)
    print(f"[{gse}] series_matrix 样本解析")
    p = os.path.join(BASE, gse, f"{gse}_series_matrix.txt.gz")
    acc = None
    chars = []  # list of (line_idx, [values per sample])
    title = None
    with gzip.open(p, 'rt', errors='ignore') as f:
        idx = 0
        for line in f:
            if line.startswith('!Sample_geo_accession'):
                acc = line.rstrip('\n').split('\t')[1:]
            elif line.startswith('!Sample_title'):
                title = line.rstrip('\n').split('\t')[1:]
            elif line.startswith('!Sample_characteristics_ch1'):
                vals = line.rstrip('\n').split('\t')[1:]
                chars.append((idx, vals))
            idx += 1
    print(f"  样本数(gsm)={len(acc) if acc else 0}")
    if acc:
        print(f"  前4 GSM={acc[:4]}  末2 GSM={acc[-2:]}")
    if title:
        print(f"  前4 title={title[:4]}")
    print(f"  !Sample_characteristics_ch1 共 {len(chars)} 行(每条=一个特征变量)")
    for li,(i,vals) in enumerate(chars):
        uniq = sorted(set(vals))
        print(f"    特征行{li}: 样本数={len(vals)}  唯一值(前12)={uniq[:12]}  {'...' if len(uniq)>12 else ''}")
    return acc, chars

# ---- GSE200744 ----
h1 = explore_count("GSE200744", "GSE200744_RNAseq_count_in_Mets.txt.gz")
acc1, chars1 = parse_series_samples("GSE200744")

# ---- GSE145412 ----
h2 = explore_count("GSE145412", "GSE145412_rnaseq.txt.gz")
acc2, chars2 = parse_series_samples("GSE145412")
