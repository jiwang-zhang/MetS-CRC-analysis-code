import os, gzip, urllib.request, time
GEO = "https://ftp.ncbi.nlm.nih.gov/geo"
UA = {"User-Agent": "Mozilla/5.0"}
BASE = r"D:/daixie-colorectal/19.databases"

def gpl_prefix(gpl):
    num = gpl[3:]
    return gpl[:3] + num[:-3] + "nnn"

def get(url, timeout=120):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()

def download(url, dest, max_retries=4):
    for i in range(max_retries):
        try:
            data = get(url)
            if len(data) < 1000:
                raise IOError(f"过小 {len(data)}B")
            with open(dest, "wb") as f:
                f.write(data)
            return len(data)
        except Exception as e:
            print(f"  重试{i+1}: {e}"); time.sleep(3)
    raise IOError(f"失败 {url}")

gpl = "GPL6947"
prefix = gpl_prefix(gpl)
cache = os.path.join(BASE, "_gpl_cache", f"{gpl}.annot.gz")
local = os.path.join(BASE, "GSE98895", f"{gpl}.annot.gz")
os.makedirs(os.path.dirname(cache), exist_ok=True)
url = f"{GEO}/platforms/{prefix}/{gpl}/annot/{gpl}.annot.gz"
print("URL:", url)
sz = download(url, cache)
print(f"  cache -> {os.path.basename(cache)} {sz/1e6:.1f}MB")
import shutil
shutil.copy(cache, local)
print(f"  copied -> GSE98895/{gpl}.annot.gz")

# 探列名
with gzip.open(cache, "rt", errors="ignore") as f:
    header = None
    for line in f:
        if line.startswith("!platform_table_begin"): continue
        if line.startswith("#") or line.startswith("!") or line.startswith("^"): continue
        if header is None:
            header = line.rstrip("\n").split("\t")
            print("NCOLS:", len(header))
            for i, c in enumerate(header):
                if c.lower() in ("id","gene symbol","gene_symbol","symbol","gene title","entrez gene id","gene id"):
                    print(f"  col {i}: {c}")
            break
