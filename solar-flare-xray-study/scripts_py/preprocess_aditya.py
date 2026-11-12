"""preprocess_aditya.py -- Aditya-L1 SoLEXS + HEL1OS light curves -> 1-minute CSVs (read by the R pipeline).

Inputs  (downloaded from PRADAN, kept outside OneDrive):
  C:\\solar-data\\aditya\\solexs\\...\\AL1_SLX_L1_YYYYMMDD_v1.x.zip   (daily; SDD2/*.lc.gz, 1 s counts)
  C:\\solar-data\\aditya\\hel1os\\...\\HLS_YYYYMMDD_HHMMSS_*.zip      (12 h; czt1/czt2 lightcurves, ~1 s)
Outputs (project folder):
  data/processed_adl1/solexs_1min.csv   time, sxr_cps
  data/processed_adl1/hel1os_1min.csv   time, czt_20_40, czt_40_60, czt_60_80, czt_80_150  (counts/s, CZT1+CZT2)
  data/processed_adl1/preprocess_qc.csv per-file QC
"""
import glob, gzip, io, os, re, sys, zipfile
from multiprocessing import Pool
import numpy as np, pandas as pd
from astropy.io import fits
import warnings; warnings.filterwarnings("ignore")

SOLEXS = r"C:\solar-data\aditya\solexs"
HEL1OS = r"C:\solar-data\aditya\hel1os"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "processed_adl1")
os.makedirs(OUT, exist_ok=True)
MIN_VALID_S = 30            # require >= 30 valid seconds in a minute


def per_minute(t_unix, y, min_valid=MIN_VALID_S):
    ok = np.isfinite(y) & np.isfinite(t_unix)
    t, y = t_unix[ok], y[ok]
    if t.size == 0:
        return pd.Series(dtype=float)
    s = pd.Series(y, index=(np.floor(t / 60).astype(np.int64)))
    g = s.groupby(level=0)
    return g.mean()[g.size() >= min_valid]


def solexs_day(f):
    day = re.search(r"_(\d{8})_v", f).group(1)
    z = zipfile.ZipFile(f)
    names = [n for n in z.namelist() if n.endswith(".lc.gz") and "SDD2" in n] or [n for n in z.namelist() if n.endswith(".lc.gz")]
    if not names:
        return day, None, "no_lc"
    h = fits.open(io.BytesIO(gzip.decompress(z.read(names[0]))))
    d = h["RATE"].data if "RATE" in [x.name for x in h] else h[1].data
    t = np.asarray(d["TIME"], float); c = np.asarray(d["COUNTS"], float)
    c[~np.isfinite(c)] = np.nan
    m = per_minute(t, c)
    return day, m, "ok:%d" % len(m)


def hel_segment(f):
    z = zipfile.ZipFile(f)
    base = os.path.basename(f)
    bands = {"czt_20_40": ["20.00KEV_TO_40.00KEV"], "czt_40_60": ["40.00KEV_TO_60.00KEV"],
             "czt_60_80": ["60.00KEV_TO_80.00KEV"], "czt_80_150": ["80.00KEV_TO_150.00KEV"]}
    per_det = {}
    for det in ("czt1", "czt2"):
        member = [n for n in z.namelist() if n.endswith("lightcurve_%s.fits" % det)]
        if not member:
            continue
        h = fits.open(io.BytesIO(z.read(member[0])))
        for bname, keys in bands.items():
            ext = [x for x in h if x.name.upper().startswith(det.upper() + "_LC_BAND") and any(k in x.name.upper() for k in keys)]
            if not ext:
                continue
            d = ext[0].data
            t = (np.asarray(d["MJD"], float) - 40587.0) * 86400.0
            r = np.asarray(d["CTR"], float)
            m = per_minute(t, r)
            per_det.setdefault(bname, []).append(m)
    out = {}
    for bname, ms in per_det.items():
        if len(ms) == 2:                       # require both detectors present, then sum rates
            out[bname] = ms[0].add(ms[1], fill_value=np.nan).dropna()
        else:
            out[bname] = ms[0]
    return base, out, "ok"


def newest_per_day(files, pat):
    best = {}
    for f in files:
        m = re.search(pat, os.path.basename(f))
        key = m.group(1); ver = tuple(int(x) for x in re.findall(r"\d+", m.group(2)))
        if key not in best or ver > best[key][0]:
            best[key] = (ver, f)
    return [v[1] for k, v in sorted(best.items())]


if __name__ == "__main__":
    qc = []
    # ---- SoLEXS ----
    sfiles = newest_per_day(glob.glob(os.path.join(SOLEXS, "**", "AL1_SLX_L1_*.zip"), recursive=True), r"_(\d{8})_(v[\d.]+)\.zip")
    print("SoLEXS days:", len(sfiles), flush=True)
    with Pool(6) as p:
        res = p.map(solexs_day, sfiles)
    parts = []
    for day, m, st in res:
        qc.append(dict(instrument="SoLEXS", file=day, status=st))
        if m is not None and len(m): parts.append(m)
    sx = pd.concat(parts).groupby(level=0).mean().sort_index()
    sx.index = pd.to_datetime(sx.index * 60, unit="s", utc=True)
    sx.rename("sxr_cps").rename_axis("time").reset_index().to_csv(os.path.join(OUT, "solexs_1min.csv"), index=False)
    print("SoLEXS minutes:", len(sx), sx.index.min(), "->", sx.index.max(), flush=True)

    # ---- HEL1OS ----
    hfiles = sorted(glob.glob(os.path.join(HEL1OS, "**", "HLS_*.zip"), recursive=True))
    # keep newest processing version per 12 h segment (V111 < V211 < V311)
    seg = {}
    for f in hfiles:
        m = re.search(r"HLS_(\d{8})_(\d{2})\d{4}_.*_V(\d+)\.zip", os.path.basename(f))
        k = (m.group(1), m.group(2)); v = int(m.group(3))
        if k not in seg or v > seg[k][0]: seg[k] = (v, f)
    hfiles = [v[1] for k, v in sorted(seg.items())]
    print("HEL1OS segments:", len(hfiles), flush=True)
    with Pool(6) as p:
        res = p.map(hel_segment, hfiles)
    cols = {}
    for base, out, st in res:
        qc.append(dict(instrument="HEL1OS", file=base, status="%s bands=%d" % (st, len(out))))
        for b, m in out.items():
            cols.setdefault(b, []).append(m)
    df = pd.DataFrame({b: pd.concat(ms).groupby(level=0).mean() for b, ms in cols.items()}).sort_index()
    df.index = pd.to_datetime(df.index * 60, unit="s", utc=True)
    df.rename_axis("time").reset_index().to_csv(os.path.join(OUT, "hel1os_1min.csv"), index=False)
    print("HEL1OS minutes:", len(df), df.index.min(), "->", df.index.max(), flush=True)
    print(df.describe().T[["count", "mean", "50%", "max"]].round(2).to_string(), flush=True)
    pd.DataFrame(qc).to_csv(os.path.join(OUT, "preprocess_qc.csv"), index=False)
    print("done", flush=True)
