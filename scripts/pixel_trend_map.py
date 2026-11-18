"""
Pixel-level NDVI trend map (Mann-Kendall + Sen's slope) for the SCCL study area.

What it does
  1. Builds the 5 km study area from QIGS files/Selected_Mines.gpkg (1, 3, 5 km zones).
  2. Reads QIGS files/NDVI_YYYY.tif for 2005-2026, cuts each to the 5 km area and
     puts every year on the same 30 m grid.
  3. For every pixel runs the Mann-Kendall test and Sen's slope on its 22-year series.
  4. Writes GeoTIFFs (slope, Z, p, trend class) and a per-zone summary CSV.

Run from the lulc_project folder:
  uv add rasterio geopandas scipy pandas numpy      (once)
  uv run python pixel_trend_map.py

Outputs go to:  Pixel trend results/
"""
from pathlib import Path
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.features import geometry_mask, rasterize
from rasterio.warp import reproject, Resampling
from scipy.stats import norm

# ---------------- settings ----------------
ROOT = Path(__file__).resolve().parent
NDVI_DIR = ROOT / "QIGS files"
MINES = NDVI_DIR / "Selected_Mines.gpkg"
OUT = ROOT / "Pixel trend results"
YEARS = list(range(2005, 2027))     # same period as the zone analysis
MIN_YEARS = 18                       # skip pixels with fewer valid years (Landsat 7 gaps, nodata)
ALPHA = 0.05
NDVI_MIN, NDVI_MAX = -1.0, 1.0       # anything outside this is treated as nodata
ROWS_PER_CHUNK = 100                 # lower this if you run out of memory
CRS = "EPSG:32644"                   # WGS 84 / UTM zone 44N
# ------------------------------------------


def build_zones():
    leases = gpd.read_file(MINES).to_crs(CRS)
    u = leases.union_all() if hasattr(leases, "union_all") else leases.unary_union
    b1, b3, b5 = u.buffer(1000), u.buffer(3000), u.buffer(5000)
    zones = gpd.GeoDataFrame(
        {"zone": ["0-1 km", "1-3 km", "3-5 km"], "code": [1, 2, 3]},
        geometry=[b1, b3.difference(b1), b5.difference(b3)], crs=CRS)
    return u, b5, zones


def reference_grid(b5, res=30.0):
    minx, miny, maxx, maxy = b5.bounds
    # snap to the Landsat 30 m grid (Landsat UTM corners sit on multiples of 15 m)
    minx = np.floor(minx / 30) * 30 + 15
    maxy = np.ceil(maxy / 30) * 30 + 15
    width = int(np.ceil((maxx - minx) / res))
    height = int(np.ceil((maxy - miny) / res))
    transform = rasterio.transform.from_origin(minx, maxy, res, res)
    return transform, width, height


def load_stack(transform, width, height, inside):
    stack = np.full((len(YEARS), height, width), np.nan, dtype=np.float32)
    for k, yr in enumerate(YEARS):
        f = NDVI_DIR / f"NDVI_{yr}.tif"
        if not f.exists():
            print(f"  missing {f.name}, year left empty")
            continue
        with rasterio.open(f) as src:
            dst = np.full((height, width), np.nan, dtype=np.float32)
            reproject(
                source=rasterio.band(src, 1), destination=dst,
                src_transform=src.transform, src_crs=src.crs, src_nodata=src.nodata,
                dst_transform=transform, dst_crs=CRS, dst_nodata=np.nan,
                resampling=Resampling.nearest)
        dst[(dst < NDVI_MIN) | (dst > NDVI_MAX)] = np.nan
        dst[~inside] = np.nan
        stack[k] = dst
        print(f"  {yr}: {np.isfinite(dst).sum():,} valid pixels")
    return stack


def mk_sen_chunk(x):
    """x: (n_years, n_pixels) with NaN for missing. Returns slope, z, p, tau, n_valid."""
    n_t, n_p = x.shape
    t = np.arange(n_t, dtype=np.float32)
    i, j = np.triu_indices(n_t, k=1)
    d = x[j] - x[i]                                  # (pairs, pixels)
    ok = np.isfinite(d)
    s = np.where(ok, np.sign(d), 0).sum(axis=0)
    n = np.isfinite(x).sum(axis=0).astype(np.float64)
    var = n * (n - 1) * (2 * n + 5) / 18.0           # NDVI is continuous, ties are negligible
    z = np.where(s > 0, (s - 1) / np.sqrt(var), np.where(s < 0, (s + 1) / np.sqrt(var), 0.0))
    p = 2 * (1 - norm.cdf(np.abs(z)))
    tau = s / (n * (n - 1) / 2)
    slopes = np.where(ok, d / (t[j] - t[i])[:, None], np.nan)
    slope = np.nanmedian(slopes, axis=0)
    bad = n < MIN_YEARS
    for a in (slope, z, p, tau):
        a[bad] = np.nan
    return slope, z, p, tau, n


def write_tif(path, arr, transform, dtype="float32", nodata=np.nan):
    with rasterio.open(path, "w", driver="GTiff", height=arr.shape[0], width=arr.shape[1],
                       count=1, dtype=dtype, crs=CRS, transform=transform,
                       nodata=nodata, compress="deflate") as dst:
        dst.write(arr.astype(dtype), 1)


def main():
    OUT.mkdir(exist_ok=True)
    print("Building zones...")
    leases, b5, zones = build_zones()
    transform, width, height = reference_grid(b5)
    inside = ~geometry_mask([b5], out_shape=(height, width), transform=transform)
    zone_r = rasterize([(g, c) for g, c in zip(zones.geometry, zones.code)],
                       out_shape=(height, width), transform=transform, fill=0, dtype="uint8")
    lease_r = rasterize([(leases, 1)], out_shape=(height, width), transform=transform,
                        fill=0, dtype="uint8")
    print(f"Grid {width} x {height}, {inside.sum():,} pixels inside 5 km")

    print("Reading NDVI rasters...")
    stack = load_stack(transform, width, height, inside)

    print("Running Mann-Kendall and Sen's slope per pixel...")
    slope = np.full((height, width), np.nan, np.float32)
    z = slope.copy(); p = slope.copy(); tau = slope.copy()
    nval = np.zeros((height, width), np.uint8)
    for r0 in range(0, height, ROWS_PER_CHUNK):
        r1 = min(r0 + ROWS_PER_CHUNK, height)
        block = stack[:, r0:r1, :].reshape(len(YEARS), -1)
        keep = inside[r0:r1, :].ravel()
        out = [np.full(keep.size, np.nan, np.float32) for _ in range(4)]
        if keep.any():
            res = mk_sen_chunk(block[:, keep])
            for o, rr in zip(out, res[:4]):
                o[keep] = rr
            nv = np.zeros(keep.size); nv[keep] = res[4]
            nval[r0:r1, :] = nv.reshape(r1 - r0, width)
        for arr, o in zip((slope, z, p, tau), out):
            arr[r0:r1, :] = o.reshape(r1 - r0, width)
        print(f"  rows {r0}-{r1} of {height}")

    # trend class: 1 sig. decrease, 2 no sig. trend, 3 sig. increase, 0 nodata
    cls = np.zeros((height, width), np.uint8)
    valid = np.isfinite(slope)
    cls[valid] = 2
    cls[valid & (p < ALPHA) & (slope < 0)] = 1
    cls[valid & (p < ALPHA) & (slope > 0)] = 3

    write_tif(OUT / "sen_slope_per_year.tif", slope, transform)
    write_tif(OUT / "sen_slope_per_decade.tif", slope * 10, transform)
    write_tif(OUT / "mk_z.tif", z, transform)
    write_tif(OUT / "mk_p.tif", p, transform)
    write_tif(OUT / "kendall_tau.tif", tau, transform)
    write_tif(OUT / "trend_class.tif", cls, transform, dtype="uint8", nodata=0)
    write_tif(OUT / "valid_years.tif", nval, transform, dtype="uint8", nodata=0)

    # per-zone summary; lease interiors reported separately from the land outside them
    rows = []
    groups = [("0-1 km (all)", zone_r == 1),
              ("0-1 km inside leases", (zone_r == 1) & (lease_r == 1)),
              ("0-1 km outside leases", (zone_r == 1) & (lease_r == 0)),
              ("1-3 km", zone_r == 2), ("3-5 km", zone_r == 3),
              ("Whole 5 km area", inside)]
    for name, m in groups:
        m = m & valid
        npx = int(m.sum())
        if npx == 0:
            continue
        rows.append({
            "Zone": name, "Pixels": npx, "Area_km2": round(npx * 0.0009, 1),
            "Sig_increase_%": round(100 * (cls[m] == 3).mean(), 1),
            "No_sig_trend_%": round(100 * (cls[m] == 2).mean(), 1),
            "Sig_decrease_%": round(100 * (cls[m] == 1).mean(), 1),
            "Median_slope_per_yr": round(float(np.nanmedian(slope[m])), 5),
            "Mean_slope_per_yr": round(float(np.nanmean(slope[m])), 5),
        })
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT / "pixel_trend_summary_by_zone.csv", index=False)
    print("\n", summary.to_string(index=False))
    print(f"\nDone. Files in: {OUT}")


if __name__ == "__main__":
    main()
