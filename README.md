# Thesis Analysis: Moisture Sources and Precipitation Extremes (July 1997 CZ Flood)

This repository contains the R scripts and outputs used for my Bachelor's thesis on identifying moisture sources and precipitation anomalies contributing to the July 1997 flood in Central Europe. Specifically, the Czech Republic was analyzed.

## Structure

- **scripts/** — R scripts for land/ocean source analysis, basin correlations, and exceedance probabilities
- **outputs/** — Figures and Excel tables from the analyses
- **data/** — Required input data (not all raw datasets were included due to size limitations)

## Scripts

- `E2P_EPs_land_ocean_country_analysis.R`: Calculates daily moisture contributions from land and ocean and by country.
- `Spearman_correlation_analysis.R`: Computes rolling Spearman correlations between basin-level evaporation and precipitation.
- `exceedance_probabilities.R`: Calculates precipitation exceedance probabilities for 14 and 31-day sums against a 1980–2020 climatology.

## Outputs

All figures and summary tables used in the thesis are saved in `/outputs`, including:
- Land vs ocean contribution time series
- Basin-level correlation plots
- Precipitation exceedance maps and tables
