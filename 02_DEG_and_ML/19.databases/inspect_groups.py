import gzip, os

BASE = r"D:/daixie-colorectal/19.databases"

def inspect(gse):
    p = os.path.join(BASE, gse, f"{gse}_series_matrix.txt.gz")
    meta = {}
    with gzip.open(p, "rt", errors="ignore") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!Sample_") or line.startswith("!series_"):
                parts = line.split("\t")
                meta[parts[0]] = parts[1:]
    ns = len(meta.get("!Sample_title", []))
    print(f"\n===== {gse}  (nsamples={ns}) =====")
    for k in ["!Sample_title", "!Sample_source_name_ch1", "!Sample_characteristics_ch1", "!Sample_description"]:
        if k in meta:
            vals = meta[k]
            # show distinct short values
            uniq = []
            for v in vals:
                v = v.strip('"')
                if len(v) > 80:
                    v = v[:80] + "..."
                if v not in uniq:
                    uniq.append(v)
            print(f"  {k}: n_unique={len(uniq)}")
            for u in uniq[:8]:
                print("     ", repr(u))

for g in ["GSE98895", "GSE21510", "GSE32323", "GSE39582"]:
    inspect(g)

# GPL570.annot header
print("\n===== GPL570.annot.gz header =====")
with gzip.open(os.path.join(BASE, "GSE98895", "GPL570.annot.gz"), "rt", errors="ignore") as f:
    for i, line in enumerate(f):
        if i < 3:
            print(line.rstrip("\n")[:300])
        else:
            break
