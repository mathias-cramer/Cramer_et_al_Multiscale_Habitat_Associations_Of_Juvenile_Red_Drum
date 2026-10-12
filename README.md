# Multiscale Habitat Associations of Juvenile Red Drum *Sciaenops ocellatus* in Tidal Rivers of Tampa Bay, Florida

## General

This repository contains the data, GIS materials, statistical analysis scripts, and results tables associated with the manuscript:

**Cramer, M. T., Radabaugh, K. R., Hearne, S. N., Adams, C. S., Stevens, P. W., Blewett, D. A., & Moyer, R. P.**  
*Multiscale Habitat Associations of Juvenile Red Drum Sciaenops ocellatus in Tidal Rivers of Tampa Bay, Florida.*

The study evaluates habitat associations of juvenile Red Drum in the Alafia, Little Manatee, and Manatee rivers of Tampa Bay, Florida. Detailed shoreline habitat assessments conducted in 2019 were integrated with fisheries-independent monitoring data collected from 2014 through 2024. Analyses focused on two juvenile size classes, 0–50 mm and 51–100 mm standard length, to evaluate how habitat associations differed across early ontogeny.

The repository is intended to provide the materials necessary to reproduce the statistical analyses and principal results reported in the manuscript. All study materials are contained within the `Materials` folder.

---

## Contents

| File | Description |
|---|---|
| `Tampa_Bay_Red_Drum_Data.csv` | Final site-level analysis dataset used throughout the study. Each row represents one of 51 shoreline assessment sites in the Alafia, Little Manatee, or Manatee rivers. The dataset contains juvenile Red Drum abundance metrics, shoreline and channel habitat characteristics, aggregated fisheries-independent monitoring environmental variables, substrate and sediment measurements, and vegetation percent-cover data. The primary manuscript analyses use the 0–50 mm and 51–100 mm standard-length abundance variables. |
| `Seine_Transect_Distance_Raw.csv` | Spatial-association table containing the 492 qualifying seine–shoreline associations within 250 m of the 51 shoreline assessment sites. The file includes seine and shoreline identifiers, separation distances, and nearest-neighbor ranks and is used to calculate median seine-to-transect distance for the spatial sensitivity analysis. |
| `GIS_PACKAGE_FILENAME` | GIS data associated with the shoreline assessment sites and spatial components of the study. |
| `Results Tables.xlsx` | Workbook containing the principal statistical results tables and supplementary tables. Sheets include model coefficient tables and Type II likelihood-ratio test results for the 0–50 mm and 51–100 mm Red Drum models (Tables 2–5), as well as complete indicator-species results for both size classes (Supplementary Tables S1–S2). |
| `Tampa Bay Red Drum, Correlation Analysis.R` | Calculates Spearman rank correlations among continuous habitat and environmental variables and juvenile Red Drum abundance, identifies strong correlations, and generates the correlation matrix presented as Figure 2. |
| `Tampa Bay Red Drum, GLM Analysis.R` | Primary modeling script. Performs data screening and collinearity assessment, Box–Cox transformation of juvenile abundance, Gaussian generalized linear modeling for both size classes, AIC-based model selection, Type II likelihood-ratio tests, Tukey-adjusted pairwise comparisons of estimated marginal means, model diagnostics, and generation of model-prediction figures. |
| `Tampa Bay Red Drum, Plant Indicator Analysis.R` | Performs permutation-based indicator-species analyses of shoreline vegetation for the 0–50 mm and 51–100 mm size classes. The script calculates indicator values, specificity, fidelity, mean vegetation cover, raw and Benjamini–Hochberg-adjusted *p*-values, and generates the indicator-species figures and results used in Supplementary Tables S1–S2. |
| `Tampa Bay Red Drum, Mangrove Analysis.R` | Conducts exploratory analyses restricted to mangrove-dominated shoreline sites. Analyses include high- versus low-abundance habitat comparisons and Spearman rank correlations between continuous habitat characteristics and juvenile Red Drum abundance. |
| `Tampa Bay Red Drum, Post Hoc Exploratory Analysis.R` | Contains secondary exploratory and post hoc analyses used to further interpret patterns identified by the primary models, including comparisons among rivers, shore types, bottom types, and embayment categories and model-adjusted comparisons of categorical predictors. |
| `Tampa Bay Red Drum, Sensitivity Analysis.R` | Evaluates whether spatial separation between fisheries-independent seine samples and shoreline assessment sites influenced model inference. The final model for each size class is refit with median seine-to-transect distance included as an additional covariate and compared with the corresponding primary model. |

---

## Reproducibility and Run Order

Analyses were conducted in **R 4.4.1**, and spatial analyses and map preparation were conducted using **ArcGIS Pro 3.1.3**.

The R scripts use relative file names and are designed to be run with the working directory set to the `Materials` folder.

For a complete analysis run within a single R session, the recommended order is:

1. `Tampa Bay Red Drum, GLM Analysis.R`
2. `Tampa Bay Red Drum, Correlation Analysis.R`
3. `Tampa Bay Red Drum, Plant Indicator Analysis.R`
4. `Tampa Bay Red Drum, Mangrove Analysis.R`
5. `Tampa Bay Red Drum, Post Hoc Exploratory Analysis.R`
6. `Tampa Bay Red Drum, Sensitivity Analysis.R`

`Tampa Bay Red Drum, GLM Analysis.R` should be run before the Plant Indicator, Mangrove, and Post Hoc Exploratory scripts because it creates shared data and model objects used by those analyses.

The Correlation Analysis script imports `Tampa_Bay_Red_Drum_Data.csv` independently and may therefore be run separately.

The Sensitivity Analysis script may also be run independently, provided that both `Tampa_Bay_Red_Drum_Data.csv` and `Seine_Transect_Distance_Raw.csv` are available in the working directory.

Required R packages are identified within the individual scripts.

---

## Data Availability, Licensing, and Use

The materials provided in this repository are publicly available for purposes including reproducibility, research, education, and further scientific analysis.

The repository contains the final analysis-ready site-level dataset and derived spatial-association data used for the analyses reported in the manuscript. It does not reproduce the complete underlying fisheries-independent monitoring database or all intermediate working files used during data preparation.

### Code

The R scripts in this repository are distributed under the **MIT License**. See [`LICENSE_CODE`](LICENSE_CODE) for details.

### Data and other research materials

The data, GIS materials, and tabular results in the `Materials` folder are distributed under the **Creative Commons Attribution 4.0 International License (CC BY 4.0)**. See [`LICENSE_DATA`](LICENSE_DATA) for details.

---

## Citation and Contact

The formal article citation and DOI will be added to this repository following publication.

Questions regarding the repository, data, or analyses may be directed to the corresponding author:

**Mathias T. Cramer**  
Florida Fish and Wildlife Conservation Commission  
Fish and Wildlife Research Institute
