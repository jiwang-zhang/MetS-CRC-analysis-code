import re, csv, os, statistics

base = r'D:\BaiduSyncdisk\daixie-colorectal\第二部分DEG\7.databases - 副本\GSE39582\GSE39582_series_matrix.txt\GSE39582_series_matrix.txt'
outdir = r'D:\BaiduSyncdisk\daixie-colorectal\第二部分DEG\26.prognostic_model'
out = os.path.join(outdir, 'GSE39582_survival.csv')

samples = None
fields = {}

with open(base, 'r', encoding='utf-8', errors='ignore') as f:
    for line in f:
        if line.startswith('!Sample_geo_accession'):
            parts = line.rstrip('\n').split('\t')
            samples = [p.strip('"') for p in parts[1:]]
        elif line.startswith('!Sample_characteristics_ch1'):
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            first_val = parts[1].strip('"')
            m = re.match(r'^([a-zA-Z_][a-zA-Z0-9_\.]*)\s*(\([^)]*\))?\s*:\s*(.*)$', first_val)
            if m:
                key = m.group(1) + (m.group(2) or '')
                vals = [p.strip('"') for p in parts[1:]]
                fields[key] = vals

def get(opts):
    for k in opts:
        if k in fields:
            return fields[k]
    return None

os_event = get(['os.event'])
os_time  = get(['os.delay(months)', 'os.delay (months)'])
rfs_event= get(['rfs.event'])
rfs_time = get(['rfs.delay', 'rfs.delay(months)', 'rfs.delay (months)'])
age      = get(['age.at.diagnosis(year)'])
sex      = get(['Sex'])
stage    = get(['tnm.stage'])
mmr      = get(['mmr.status'])

def split_val(v):
    if v is None:
        return ''
    return v.split(':')[-1].strip()

def to_num(v):
    try:
        return float(v)
    except:
        return ''

with open(out, 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['Sample','os_event','os_time_months','rfs_event','rfs_time_months',
                'age_year','sex','tnm_stage','mmr_status'])
    for i, s in enumerate(samples):
        oe = split_val(os_event[i]) if os_event else ''
        ot = to_num(split_val(os_time[i])) if os_time else ''
        re_= split_val(rfs_event[i]) if rfs_event else ''
        rt = to_num(split_val(rfs_time[i])) if rfs_time else ''
        a  = split_val(age[i]) if age else ''
        sx = split_val(sex[i]) if sex else ''
        st = split_val(stage[i]) if stage else ''
        mm = split_val(mmr[i]) if mmr else ''
        w.writerow([s, oe, ot, re_, rt, a, sx, st, mm])

print('Saved:', out)
print('Total samples:', len(samples))

os_t  = [float(split_val(os_time[i])) for i in range(len(samples)) if os_time and split_val(os_time[i]) not in ('', 'NA', 'N/A')]
os_e  = [int(split_val(os_event[i])) for i in range(len(samples)) if os_event and split_val(os_event[i]) in ('0','1')]
rfs_t = [float(split_val(rfs_time[i])) for i in range(len(samples)) if rfs_time and split_val(rfs_time[i]) not in ('', 'NA', 'N/A')]
rfs_e = [int(split_val(rfs_event[i])) for i in range(len(samples)) if rfs_event and split_val(rfs_event[i]) in ('0','1')]

print(f'OS : n={len(os_t)}, events={sum(os_e)}, median={statistics.median(os_t):.1f}mo, max={max(os_t):.1f}mo')
print(f'RFS: n={len(rfs_t)}, events={sum(rfs_e)}, median={statistics.median(rfs_t):.1f}mo, max={max(rfs_t):.1f}mo')
print('First 3 rows:')
with open(out) as f:
    for j, line in enumerate(f):
        if j < 3:
            print('  ', line.strip()[:130])
