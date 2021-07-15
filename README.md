# pretermBirth_metabolomics

Data and code for the publication Link "Maternal lipids are involved in the pathogenesis of preterm labor" authored by Yile Chen, Bing He, and Yu Liu et al.

### Data
* quant.metabolitemeasurements.rds: quantile normalized metabolomics data
* clinic.rds: clinical factors

### Usage
There are six trunks in this repositories:
* pretermmetabolite_sov.R: Data preprocessing (imputation, quantile normalization and log-transformation) and SOV analysis. 
* Diff_analysis.R: differential lipid identification and marker selection
* WGCNA.R: Weighted Gene Co-Expression Network Analysis
* ML_main.R: machine learning model construction (source machine_lerning11.R
first)
* causality analysis.R: causality analysis
