# CPEC satellite analysis

Local-effects analysis of CPEC projects using Google Earth Engine exports
(nightlights around projects and roads, district built-up area, Sentinel-5P NO₂ around power plants).

**Start here:** [`docs/AUDIT.md`](docs/AUDIT.md) — data audit findings and what they change in the research plan.

## Layout

```
data/raw/        GEE CSV exports (not committed — place the five files here)
                   1_ntl_district_annual.csv  2_ntl_rings_projects.csv  3_ntl_rings_roads.csv
                   4_builtup_district.csv     6_no2_monthly.csv
data/crosswalk/  templates to fill: project UID → project/dates, 119 districts → 86 LFS units
src/01_audit.py            data dictionary, coverage, sensor break, duplicate-site clusters, raw plots
src/02_no2_diagnostics.py  core − control NO₂ by site, pre-period inventory, deseasonalised contrasts
output/          tables and figures written by the scripts
```

## Run

```
pip install -r requirements.txt
python src/01_audit.py
python src/02_no2_diagnostics.py
```
