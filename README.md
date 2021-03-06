# Aging_2021
Source code accompanying the manuscript entitled "Integration of segmented regression analysis with weighted gene correlation network analysis identifies genes whose expression is remodeled throughout physiological aging in mouse tissues", published in *Aging (Albany NY)* (https://doi.org/10.18632/aging.203379).

## Author information
Margarida Ferreira<sup>1</sup>, Stephany Francisco<sup>1</sup>, Ana R. Soares<sup>1</sup>, Ana Nobre<sup>1</sup>, Miguel Pinheiro<sup>1</sup>, Andreia Reis<sup>1</sup>, Sonya Neto<sup>1</sup>, Ana João Rodrigues<sup>2,3</sup>, Nuno Sousa<sup>2,3</sup>, Gabriela Moura<sup>1</sup> and Manuel A.S. Santos<sup>1#</sup>

<sup>1</sup> iBiMED - Institute of Biomedicine, Department of Medical Sciences, University of Aveiro, 3810-193 Aveiro, Portugal

<sup>2</sup> Life and Health Sciences Research Institute (ICVS), School of Medicine, University of Minho, 4710-057 Braga, Portugal

<sup>3</sup> ICVS/3B’s–PT Government Associate Laboratory, Braga/Guimarães, Portugal

<sup>#</sup>corresponding author: Manuel A.S. Santos (msantos@ua.pt)

## Folder and file description

1. Normalization: Data pre-processing and normalization
   1. Normalization.md: markdown documentation file
   2. annotation.txt: ensembl biotype annotation; downloaded in March 25th 2021

Note: Because we are not the owners of the dataset, we do not provide here the gene expression data nor the metadata. For retrieving these files we advise you to download them directly from the GSE132040 entry of the Gene Expression Omnibus database (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132040).

2. Regression: Segmented regression analysis
   1. segmented_regression.md: markdown documentation file
   2. annotation.txt: ensembl biotype annotation; downloaded in March 25th 2021 (resulting from the pre-processing and normalization steps)
   3. norm_XXXX_counts.txt: normalized count data for each tissue (resulting from the pre-processing and normalization steps)
   4. coldata_XXXX.txt: sample information for each tissue (resulting from the pre-processing and normalization steps)
   5. res_XXXX_cov.rds: R objects corresponding to the results of the trendy function for each tissue
   6. res_XXXX_r2.rds: R objects corresponding to the choice of adjusted R2 threshold for each tissue
   7. res.top_XXXX_cov.rds: R objects corresponding to the results of the topTrendy function (top dynamic genes) for each tissue

3. WGCNA: network construction
   1. WGCNA.md: markdown documentation file
   2. vst_XXXX_counts.txt: vst-transformed count data for each tissue (resulting from the pre-processing and normalization steps)
   3. coldata_XXXX.txt: sample information for each tissue (resulting from the pre-processing and normalization steps)
  
4. WGCNA_2: module-trait associations
   1. WGCNA_2.md: markdown documentation file
   2. XXXX_dataInput.RData: RData files containing the expression and trait data for each tissue (resulting from the network construction)
   3. XXXX_networkConstruction-auto.RData: RData files containing the network data for each tissue (resulting from the network construction)
   4. XXXX_table_sig.txt: top dynamic genes for each tissue (resulting from the segmented regression analysis)

5. Functional_analysis: gene ontology over-representation analysis
   1. functional_analysis.md: markdown documentation file
   2. vst_XXXX_counts.txt: vst-transformed count data for each tissue (resulting from the pre-processing and normalization steps)
   3. upsetList_XXXX_trendyYYYY.txt: gene lists resulting from the overlap between trendy and hub genes for each significant module of each tissue (resulting from the module-trait association analysis)
