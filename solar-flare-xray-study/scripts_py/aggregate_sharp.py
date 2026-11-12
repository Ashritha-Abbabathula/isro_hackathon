"""aggregate_sharp.py -- full-disk aggregates of SDO/HMI SHARP keywords (pure python, no pandas).
Reads data/processed_adl1/sharp_rows.csv (per region, per 12 min). Fetches any missing day via pull_sharp.fetch_day, then
writes sharp_fulldisk_12min.csv. Good-quality regions only (QUALITY == 0x00000000) within +-60 deg of the central meridian
(CEA maps are unreliable near the limb)."""
import csv, os, sys, datetime as dt
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pull_sharp import fetch_day, KEYS, OUT

rows_path = os.path.join(OUT, "sharp_rows.csv")
rows = list(csv.DictReader(open(rows_path)))
have = {r["T_REC"][:10] for r in rows}
day = dt.date(2024, 7, 1); missing = []
while day <= dt.date(2024, 12, 31):
    if day.strftime("%Y.%m.%d") not in have: missing.append(day)
    day += dt.timedelta(days=1)
print("days missing from raw rows:", [str(d) for d in missing], flush=True)
for d in missing:
    got = fetch_day(d)
    if got:
        rows += got; print("  fetched", d, len(got), flush=True)
rows.sort(key=lambda r: (r["T_REC"], r["HARPNUM"]))
with open(rows_path, "w", newline="") as f:
    w = csv.writer(f); w.writerow(KEYS)
    for r in rows: w.writerow([r[k] for k in KEYS])

def num(x):
    try:
        v = float(x); return v if v == v else None
    except Exception:
        return None

agg = defaultdict(lambda: defaultdict(list)); used = 0
for r in rows:
    try: q = int(r["QUALITY"], 16)
    except Exception: continue
    lon = num(r["LON_FWT"])
    if q != 0 or lon is None or abs(lon) > 60: continue
    used += 1; a = agg[r["T_REC"]]
    for k in ["USFLUX", "TOTUSJH", "TOTPOT", "SAVNCPP", "ABSNJZH", "R_VALUE", "SHRGT45", "AREA_ACR"]:
        v = num(r[k])
        if v is not None: a[k].append(v)
out = os.path.join(OUT, "sharp_fulldisk_12min.csv")
with open(out, "w", newline="") as f:
    w = csv.writer(f); w.writerow(["T_REC", "n_regions", "sum_usflux", "sum_totusjh", "sum_totpot", "sum_savncpp", "sum_absnjzh", "max_r_value", "mean_shrgt45", "sum_area"])
    for t in sorted(agg):
        a = agg[t]
        w.writerow([t, len(a["USFLUX"]), sum(a["USFLUX"]), sum(a["TOTUSJH"]), sum(a["TOTPOT"]), sum(a["SAVNCPP"]), sum(a["ABSNJZH"]),
                    max(a["R_VALUE"]) if a["R_VALUE"] else "", (sum(a["SHRGT45"]) / len(a["SHRGT45"])) if a["SHRGT45"] else "", sum(a["AREA_ACR"])])
print("regions-rows used:", used, "| 12-min steps written:", len(agg), flush=True)
