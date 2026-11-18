# Vegetation Dynamics under Mining and Monsoon Variability

Temporal and spatial study of the Singareni Coalfields, Telangana (M.Tech Dissertation, Phase I).

This project tracks dry-season vegetation around nine Singareni Collieries Company Limited (SCCL) mining leases between Ramagundam and Mandamarri from 2005 to 2026. Mean NDVI from one January or February Landsat 5/7/8/9 scene per year is compared across three distance zones (0–1, 1–3 and 3–5 km from the leases). It is then tested against IMD rainfall, ERA5 temperature and SCCL coal production.

Main results:

- NDVI increased in all three zones (Mann–Kendall τ = 0.57–0.65, p ≤ 0.0003; Sen's slope 0.0024–0.0026 per year).
- The 0–1 km zone, which contains the leases, was the least green in 19 of 22 years.
- Rainfall and temperature showed no significant trend.
- Raw correlations between NDVI and opencast production (r ≈ 0.70–0.76) disappear after removing the shared time trend from both series (|r| ≤ 0.28, not significant).

## Datasets are not included

The datasets are **not available in this repository** because of their size. This covers the Landsat scenes and NDVI rasters, the IMD NetCDF files, the ERA5 GRIB file, the SCCL annual reports, the lease boundary layers, and the derived CSV and figure outputs. The code expects them in the original local project folder.

Data sources:

| Data | Source |
|---|---|
| Landsat Collection 2 Level-2 surface reflectance, path 143 row 47 | USGS EarthExplorer |
| Daily gridded rainfall, 0.25° | India Meteorological Department (IMD) |
| ERA5 2 m air temperature | Copernicus Climate Data Store (ECMWF) |
| Coal production, FY 2004–05 to 2024–25 | SCCL annual reports |
| Mining lease boundaries | Telangana State Remote Sensing Applications Centre (TGRAC) |

## Repository layout

```
notebooks/
  01_imd_rainfall_extraction.ipynb          IMD NetCDF -> annual and monsoon rainfall
  02_era5_temperature_extraction.ipynb      ERA5 GRIB -> annual and Jan-Feb temperature
  03_rainfall_temperature_anomalies.ipynb   Rainfall and temperature anomalies
  04_ndvi_trend_mann_kendall_sen.ipynb      Zone NDVI trends (Mann-Kendall, Sen's slope)
  05_ndvi_anomalies.ipynb                   NDVI anomalies (2005-2025 reference)
  06_ndvi_climate_mining_correlation.ipynb  NDVI vs climate and production correlations
  07_ndvi_spatial_maps.ipynb                Spatial NDVI map panels
scripts/
  trend_and_detrended_correlation.py        Trend tests and raw vs detrended correlations
  pixel_trend_map.py                        Pixel-level Mann-Kendall / Sen's slope maps
docs/
  workflow_flowchart.drawio                 Project workflow (open in draw.io)
src/lulc_project/                           Package stub created by uv
pyproject.toml, uv.lock, .python-version    Python 3.13 project files (uv)
requirements.txt                            Packages used by the notebooks and scripts
```

NDVI rasters and zonal statistics were produced in QGIS, so that step has no code here.

## Setup

With [uv](https://docs.astral.sh/uv/):

```
uv sync
uv pip install -r requirements.txt
```

Or with pip:

```
pip install -r requirements.txt
```

The notebooks use absolute Windows paths to the original data folder (`C:\Users\RUTHVIK REDDY\lulc_project\...`). Change those paths to your own data location before running them.
