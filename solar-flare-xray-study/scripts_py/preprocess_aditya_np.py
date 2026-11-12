"""preprocess_aditya_np.py -- numpy-only version (no pandas) of the Aditya-L1 preprocessing, plus sub-minute HXR features.
Inputs : C:\\solar-data\\aditya\\solexs\\**\\AL1_SLX_L1_*.zip   and   C:\\solar-data\\aditya\\hel1os\\**\\HLS_*.zip
Outputs: <OUT>/solexs_1min.csv  (time, sxr_cps)
         <OUT>/hel1os_1min.csv  (time, czt_20_40, czt_40_60, czt_60_80, czt_80_150,
                                  lo_max8s, lo_cv8s, lo_nburst, hi_max8s, hi_cv8s, hi_nburst)
Minute definition is identical to the original pipeline: per detector, per-minute mean of the 1 s rates (>= 30 valid seconds),
then the two CZT detectors are summed. Sub-minute features use 8-second blocks (the true readout cadence) of the CZT1+CZT2 sum (both detectors required).
Usage: python preprocess_aditya_np.py <OUT_DIR> [date_from YYYYMMDD] [date_to YYYYMMDD]
"""
import glob, gzip, io, os, re, sys, zipfile, datetime as dt
from multiprocessing import Pool
import numpy as np
from astropy.io import fits
import warnings; warnings.filterwarnings("ignore")

SOLEXS = r"C:\solar-data\aditya\solexs"; HEL1OS = r"C:\solar-data\aditya\hel1os"; MIN_VALID = 30
BANDS = {"czt_20_40": "20.00KEV_TO_40.00KEV", "czt_40_60": "40.00KEV_TO_60.00KEV", "czt_60_80": "60.00KEV_TO_80.00KEV", "czt_80_150": "80.00KEV_TO_150.00KEV"}

def minute_means(t, y):
    """t unix seconds (float), y values -> dict minute_index -> mean, only minutes with >= MIN_VALID finite samples."""
    ok = np.isfinite(t) & np.isfinite(y); t, y = t[ok], y[ok]
    if t.size == 0: return {}
    m = np.floor(t / 60).astype(np.int64); u, inv = np.unique(m, return_inverse=True)
    cnt = np.bincount(inv); sm = np.bincount(inv, weights=y)
    return {int(a): float(s / c) for a, c, s in zip(u, cnt, sm) if c >= MIN_VALID}

def solexs_day(f):
    z = zipfile.ZipFile(f)
    names = [n for n in z.namelist() if n.endswith(".lc.gz") and "SDD2" in n] or [n for n in z.namelist() if n.endswith(".lc.gz")]
    if not names: return os.path.basename(f), {}
    h = fits.open(io.BytesIO(gzip.decompress(z.read(names[0]))))
    d = h["RATE"].data if "RATE" in [x.name for x in h] else h[1].data
    return os.path.basename(f), minute_means(np.asarray(d["TIME"], float), np.asarray(d["COUNTS"], float))

BLOCK = 8   # HEL1OS Level-1 'CTR' is stored per second but is a raw count total delivered about every 7-9 s (83% of samples are zero,
            # CTR == STAT_ERR**2). The real cadence is ~8 s, so sub-minute features are built on 8 s blocks, not on 1 s samples.

def block_counts(t, y):
    """Sum counts in 8 s blocks (aligned to unix time); a block needs >= 6 of its 8 one-second samples."""
    ok = np.isfinite(t) & np.isfinite(y); t, y = t[ok], y[ok]
    b = np.floor(t / BLOCK).astype(np.int64); u, inv = np.unique(b, return_inverse=True)
    n = np.bincount(inv); sm = np.bincount(inv, weights=y); keep = n >= 6
    return u[keep], sm[keep]

def per_block(det_arrays):
    """Two detectors -> (block index, total counts per second) where both detectors have a valid block."""
    (b1, c1), (b2, c2) = [block_counts(t, r) for t, r in det_arrays]
    common = np.intersect1d(b1, b2)
    return common, (c1[np.searchsorted(b1, common)] + c2[np.searchsorted(b2, common)]) / BLOCK

def sub_minute(blk, rate):
    """Per minute from 8 s blocks: max block rate, coefficient of variation, and number of blocks above median + 3*Poisson sd."""
    out = {}
    m = (blk * BLOCK) // 60; u, start = np.unique(m, return_index=True); end = np.append(start[1:], len(m))
    for a, s0, s1 in zip(u, start, end):
        r = rate[s0:s1]
        if r.size < 5: continue
        med = float(np.median(r)); mean = float(r.mean())
        out[int(a)] = (float(r.max()), float(r.std() / mean) if mean > 0 else np.nan, int(np.sum(r > med + 3.0 * np.sqrt(max(med, 0.125) / BLOCK))))
    return out

def hel_segment(f):
    z = zipfile.ZipFile(f); res = {}
    arrs = {}                                          # (det, band) -> (t, r)
    for det in ("czt1", "czt2"):
        member = [n for n in z.namelist() if n.endswith("lightcurve_%s.fits" % det)]
        if not member: return os.path.basename(f), {}, {}
        h = fits.open(io.BytesIO(z.read(member[0])))
        for b, key in BANDS.items():
            ext = [x for x in h if x.name.upper().startswith(det.upper() + "_LC_BAND") and key in x.name.upper()]
            if not ext: continue
            d = ext[0].data
            arrs[(det, b)] = ((np.asarray(d["MJD"], float) - 40587.0) * 86400.0, np.asarray(d["CTR"], float))
    mm = {}
    for b in BANDS:
        if ("czt1", b) in arrs and ("czt2", b) in arrs:
            a1 = minute_means(*arrs[("czt1", b)]); a2 = minute_means(*arrs[("czt2", b)])
            mm[b] = {k: a1[k] + a2[k] for k in a1 if k in a2}
    sm = {}
    if all(("czt1", b) in arrs and ("czt2", b) in arrs for b in BANDS):
        sl, rl = per_block([arrs[("czt1", "czt_20_40")], arrs[("czt2", "czt_20_40")]])
        s6, r6 = per_block([arrs[("czt1", "czt_60_80")], arrs[("czt2", "czt_60_80")]])
        s8, r8 = per_block([arrs[("czt1", "czt_80_150")], arrs[("czt2", "czt_80_150")]])
        common = np.intersect1d(s6, s8); hi_r = r6[np.searchsorted(s6, common)] + r8[np.searchsorted(s8, common)]
        lo = sub_minute(sl, rl); hi = sub_minute(common, hi_r)
        for k in set(lo) & set(hi): sm[k] = lo[k] + hi[k]
    return os.path.basename(f), mm, sm

def newest(files, pat):
    best = {}
    for f in files:
        m = re.search(pat, os.path.basename(f)); key = m.group(1); ver = tuple(int(x) for x in re.findall(r"\d+", m.group(2)))
        if key not in best or ver > best[key][0]: best[key] = (ver, f)
    return [v[1] for k, v in sorted(best.items())]

def iso(minute_idx): return dt.datetime.fromtimestamp(minute_idx * 60, dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

if __name__ == "__main__":
    OUT = sys.argv[1]; os.makedirs(OUT, exist_ok=True)
    d0 = sys.argv[2] if len(sys.argv) > 2 else "00000000"; d1 = sys.argv[3] if len(sys.argv) > 3 else "99999999"
    sf = [f for f in newest(glob.glob(os.path.join(SOLEXS, "**", "AL1_SLX_L1_*.zip"), recursive=True), r"_(\d{8})_(v[\d.]+)\.zip") if d0 <= re.search(r"_(\d{8})_", f).group(1) <= d1]
    print("SoLEXS days:", len(sf), flush=True)
    with Pool(6) as p: res = p.map(solexs_day, sf)
    S = {}
    for _, m in res:
        for k, v in m.items(): S.setdefault(k, []).append(v)
    with open(os.path.join(OUT, "solexs_1min.csv"), "w") as f:
        f.write("time,sxr_cps\n")
        for k in sorted(S): f.write("%s,%.6g\n" % (iso(k), float(np.mean(S[k]))))
    print("SoLEXS minutes:", len(S), flush=True)
    seg = {}
    for f in glob.glob(os.path.join(HEL1OS, "**", "HLS_*.zip"), recursive=True):
        m = re.search(r"HLS_(\d{8})_(\d{2})\d{4}_.*_V(\d+)\.zip", os.path.basename(f))
        if not (d0 <= m.group(1) <= d1): continue
        k = (m.group(1), m.group(2)); v = int(m.group(3))
        if k not in seg or v > seg[k][0]: seg[k] = (v, f)
    hf = [v[1] for k, v in sorted(seg.items())]; print("HEL1OS segments:", len(hf), flush=True)
    with Pool(6) as p: res = p.map(hel_segment, hf)
    H = {}; X = {}
    for _, mm, sm in res:
        for b, d in mm.items():
            for k, v in d.items(): H.setdefault(k, {}).setdefault(b, []).append(v)
        for k, v in sm.items(): X.setdefault(k, []).append(v)
    cols = list(BANDS) + ["lo_max8s", "lo_cv8s", "lo_nburst", "hi_max8s", "hi_cv8s", "hi_nburst"]
    with open(os.path.join(OUT, "hel1os_1min.csv"), "w") as f:
        f.write("time," + ",".join(cols) + "\n")
        for k in sorted(H):
            row = [("%.6g" % np.mean(H[k][b])) if b in H[k] else "" for b in BANDS]
            x = X.get(k); xs = np.mean(np.array(x, float), axis=0) if x else [np.nan] * 6
            row += [("%.6g" % v) if np.isfinite(v) else "" for v in xs]
            f.write("%s,%s\n" % (iso(k), ",".join(row)))
    print("HEL1OS minutes:", len(H), "| with sub-minute features:", len(X), flush=True); print("done", flush=True)
