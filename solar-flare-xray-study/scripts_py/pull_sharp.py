"""pull_sharp.py -- SDO/HMI SHARP (CEA) keywords from JSOC for a date range (pure python: no pandas needed).
Query: hmi.sharp_cea_720s, 12-minute cadence, one request per day. Writes per-region rows and full-disk aggregates.
Usage: python pull_sharp.py 2024-07-01 2024-12-31
Outputs (project data/processed_adl1/): sharp_rows.csv (per region, per 12 min) and sharp_fulldisk_12min.csv (aggregates)
"""
import csv, json, os, sys, time, datetime as dt, urllib.parse, urllib.request
from collections import defaultdict

KEYS = ["HARPNUM", "T_REC", "NOAA_AR", "QUALITY", "LON_FWT", "LAT_FWT", "USFLUX", "TOTUSJH", "TOTUSJZ", "TOTPOT", "SAVNCPP", "ABSNJZH", "R_VALUE", "SHRGT45", "MEANSHR", "AREA_ACR"]
BASE = "http://jsoc.stanford.edu/cgi-bin/ajax/jsoc_info"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "processed_adl1")
os.makedirs(OUT, exist_ok=True)

def fetch_day(day):
    ds = "hmi.sharp_cea_720s[][%s_00:00:00_TAI/1d@12m]" % day.strftime("%Y.%m.%d")
    url = BASE + "?" + urllib.parse.urlencode({"ds": ds, "op": "rs_list", "key": ",".join(KEYS)})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                d = json.loads(r.read().decode("utf8"))
            if d.get("status", 0) not in (0, None) and "keywords" not in d:
                return []
            kw = {k["name"]: k["values"] for k in d.get("keywords", [])}
            n = len(kw.get("T_REC", []))
            return [{k: kw[k][i] for k in KEYS} for i in range(n)]
        except Exception as e:
            print("  retry", day, attempt + 1, str(e)[:80], flush=True); time.sleep(5 * (attempt + 1))
    print("  FAILED", day, flush=True); return None

def num(x):
    try:
        v = float(x); return v if v == v else None
    except Exception:
        return None

if __name__ == "__main__":
    d0 = dt.date.fromisoformat(sys.argv[1]); d1 = dt.date.fromisoformat(sys.argv[2])
    rows_path = os.path.join(OUT, "sharp_rows.csv"); agg = defaultdict(lambda: defaultdict(list)); failed = []
    with open(rows_path, "w", newline="") as f:
        w = csv.writer(f); w.writerow(KEYS)
        day = d0
        while day <= d1:
            rows = fetch_day(day)
            if rows is None: failed.append(str(day))
            else:
                for r in rows:
                    w.writerow([r[k] for k in KEYS])
                    # keep good-quality regions within +-60 deg of the central meridian (CEA maps are unreliable near the limb)
                    q = num(r["QUALITY"]); lon = num(r["LON_FWT"])
                    if q is not None and int(q) == 0 and lon is not None and abs(lon) <= 60:
                        a = agg[r["T_REC"]]
                        for k in ["USFLUX", "TOTUSJH", "TOTPOT", "SAVNCPP", "ABSNJZH", "R_VALUE", "SHRGT45", "AREA_ACR"]:
                            v = num(r[k])
                            if v is not None: a[k].append(v)
            print(day, "rows:", len(rows) if rows is not None else "FAILED", flush=True)
            day += dt.timedelta(days=1); time.sleep(0.5)
    out = os.path.join(OUT, "sharp_fulldisk_12min.csv")
    with open(out, "w", newline="") as f:
        w = csv.writer(f); w.writerow(["T_REC", "n_regions", "sum_usflux", "sum_totusjh", "sum_totpot", "sum_savncpp", "sum_absnjzh", "max_r_value", "mean_shrgt45", "sum_area"])
        for t in sorted(agg):
            a = agg[t]
            w.writerow([t, len(a["USFLUX"]), sum(a["USFLUX"]), sum(a["TOTUSJH"]), sum(a["TOTPOT"]), sum(a["SAVNCPP"]), sum(a["ABSNJZH"]),
                        max(a["R_VALUE"]) if a["R_VALUE"] else "", (sum(a["SHRGT45"]) / len(a["SHRGT45"])) if a["SHRGT45"] else "", sum(a["AREA_ACR"])])
    print("done; failed days:", failed, flush=True)
