"""Zone-level Mann-Kendall / Sen's slope and raw vs detrended correlations.

Needs two CSV files that are not in this repo (see README):
  ndvi_zone_means_2005_2026.csv   columns: Year, NDVI_0_1km, NDVI_1_3km, NDVI_3_5km, NDVI_overall
  NDVI_Climate_SCCL_Merged.csv    output of notebooks/06_ndvi_climate_mining_correlation.ipynb
Edit the paths below to point at them."""
import numpy as np
import pandas as pd
import pymannkendall as mk
from scipy.stats import pearsonr

ndvi = pd.read_csv("ndvi_zone_means_2005_2026.csv")
for col in ["NDVI_0_1km", "NDVI_1_3km", "NDVI_3_5km", "NDVI_overall"]:
    r = mk.original_test(ndvi[col])
    pw = mk.pre_whitening_modification_test(ndvi[col])
    print(f"{col}: tau={r.Tau:.3f} p={r.p:.5f} sen={r.slope:.6f} p_prewhitened={pw.p:.3f}")

# NDVI anomalies paired with climate and production (21 years, 2005-2025)
m = pd.read_csv("NDVI_Climate_SCCL_Merged.csv")      # from the project folder
yr = m["Year"].to_numpy(float)

def detrend(v):
    b, a = np.polyfit(yr, v, 1)          # straight line in year
    return v - (a + b * yr)              # what is left after removing it

for zone in ["0–1 km zone anomaly", "1–3 km zone anomaly", "3–5 km zone anomaly", "Overall anomaly"]:
    for drv in ["Previous_Monsoon_Rainfall_Anomaly_mm", "Jan_Feb_Temperature_Anomaly_C",
                "Open_Cast_Production", "Underground_Production"]:
        r_raw, p_raw = pearsonr(m[zone], m[drv])
        r_det, p_det = pearsonr(detrend(m[zone].to_numpy()), detrend(m[drv].to_numpy()))
        print(f"{zone:22s} {drv:38s} raw r={r_raw:+.2f} (p={p_raw:.3f})  detrended r={r_det:+.2f} (p={p_det:.3f})")
