import os, re, shutil, gzip, urllib.request, time

BASE = r"D:\daixie-colorectal\19.databases"
GEO = "https://ftp.ncbi.nlm.nih.gov/geo"

# (GSE, 预期GPL或None, 可选RNA-seq计数补充文件名)
datasets = [
    ("GSE98895",  "GPL570",    None),
    ("GSE200744", "GPL28038",  "GSE200744_RNAseq_count_in_Mets.txt.gz"),
    ("GSE145412", None,        None),   # 敏感性；计数矩阵=已下的 rnaseq.txt.gz
    ("GSE21510",  "GPL570",    None),
    ("GSE32323",  "GPL570",    None),
    ("GSE39582",  "GPL570",    None),
    ("GSE10950",  "GPL6104",   None),
    ("GSE21815",  "GPL6480",   None),
]

UA = {"User-Agent": "Mozilla/5.0 (dataset download)"}
gpl_cache = {}


def gse_prefix(gse):
    num = gse[3:]
    return gse[:3] + num[:-3] + "nnn"


def gpl_prefix(gpl):
    num = gpl[3:]
    return gpl[:3] + num[:-3] + "nnn"


def get(url, timeout=90):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def download(url, dest, max_retries=4, timeout=120):
    for i in range(max_retries):
        try:
            data = get(url, timeout)
            if len(data) < 1000:
                raise IOError(f"过小({len(data)}B)，疑似错误页")
            with open(dest, "wb") as f:
                f.write(data)
            return len(data)
        except Exception as e:
            print(f"  [重试 {i+1}/{max_retries}] {e}")
            time.sleep(3)
    raise IOError(f"下载失败: {url}")


def parse_platform(matrix_path):
    with gzip.open(matrix_path, "rt", errors="ignore") as f:
        for _ in range(120):
            line = f.readline()
            # 兼容两种写法: !Platform_ID 与 !Series_platform_id
            if line.startswith("!Platform_ID") or line.startswith("!Series_platform_id"):
                m = re.search(r"GPL\d+", line)
                if m:
                    return m.group(0)
    raise IOError("无法解析 Platform_ID")


def ensure_gpl(gpl):
    if gpl in gpl_cache:
        return gpl_cache[gpl]
    prefix = gpl_prefix(gpl)
    cache_dir = os.path.join(BASE, "_gpl_cache")
    os.makedirs(cache_dir, exist_ok=True)
    # 先试 annot.gz，失败回退 soft family
    candidates = [
        (f"{GEO}/platforms/{prefix}/{gpl}/annot/{gpl}.annot.gz",
         os.path.join(cache_dir, f"{gpl}.annot.gz")),
        (f"{GEO}/platforms/{prefix}/{gpl}/soft/{gpl}_family.soft.gz",
         os.path.join(cache_dir, f"{gpl}.soft.gz")),
    ]
    for url, dest in candidates:
        if os.path.exists(dest) and os.path.getsize(dest) > 1000:
            gpl_cache[gpl] = dest
            return dest
        try:
            print(f"[GPL] 尝试 {gpl}: {url.split('/')[-1]}")
            sz = download(url, dest)
            print(f"      -> {os.path.basename(dest)} ({sz/1e6:.1f} MB)")
            gpl_cache[gpl] = dest
            return dest
        except Exception as e:
            print(f"      失败: {e}")
    print(f"  [警告] {gpl} 注释文件两种路径均无法下载，跳过 GPL（RNA-seq 计数矩阵多自带基因名，影响有限）")
    return None


def main():
    os.makedirs(BASE, exist_ok=True)
    summary = []
    for gse, expect_gpl, suppl in datasets:
        folder = os.path.join(BASE, gse)
        os.makedirs(folder, exist_ok=True)
        print(f"\n=== {gse} ===")
        sm_name = f"{gse}_series_matrix.txt.gz"
        sm_path = os.path.join(folder, sm_name)
        if os.path.exists(sm_path) and os.path.getsize(sm_path) > 1000:
            print(f"  series_matrix 已存在 ({os.path.getsize(sm_path)/1e6:.1f} MB)")
        else:
            url = f"{GEO}/series/{gse_prefix(gse)}/{gse}/matrix/{sm_name}"
            print(f"  下载 series_matrix")
            sz = download(url, sm_path)
            print(f"      -> {sm_name} ({sz/1e6:.1f} MB)")
        gpl = expect_gpl or parse_platform(sm_path)
        print(f"  Platform_ID = {gpl}")
        gpl_src = ensure_gpl(gpl)
        if gpl_src:
            gpl_dst = os.path.join(folder, os.path.basename(gpl_src))
            if not os.path.exists(gpl_dst) or os.path.getsize(gpl_dst) < 1000:
                shutil.copy(gpl_src, gpl_dst)
            print(f"  GPL -> {os.path.basename(gpl_dst)} ({os.path.getsize(gpl_dst)/1e6:.1f} MB)")
        # RNA-seq 计数矩阵
        if suppl:
            cnt_dst = os.path.join(folder, suppl)
            if os.path.exists(cnt_dst) and os.path.getsize(cnt_dst) > 1000:
                print(f"  计数矩阵已存在 ({os.path.getsize(cnt_dst)/1e6:.1f} MB)")
            else:
                url = f"{GEO}/series/{gse_prefix(gse)}/{gse}/suppl/{suppl}"
                print(f"  下载计数矩阵: {suppl}")
                sz = download(url, cnt_dst)
                print(f"      -> {suppl} ({sz/1e6:.1f} MB)")
        summary.append(gse)

    # 整理已有的 GSE145412 计数矩阵
    loose = os.path.join(BASE, "GSE145412_rnaseq.txt.gz")
    if os.path.exists(loose):
        dst = os.path.join(BASE, "GSE145412", "GSE145412_rnaseq.txt.gz")
        shutil.move(loose, dst)
        print(f"\n[整理] 移动 GSE145412_rnaseq.txt.gz -> GSE145412/")

    print("\n========== 完成 ==========")
    print("已处理:", ", ".join(summary))


if __name__ == "__main__":
    main()
