# CNS_RF_Class
This repository contains part of the files used for the analysis presented in the master's thesis: "Effects of variable down-sampling and imputation of missing values from Nanopore sequencing on DNA methylation-based classification of tumors from the central nervous system" presented by Camilo Jose Hernandez-Toro for the degree of M.Sc. in Bioinformatics from Freie Universität Berlin.

## Central Nevous System tumor Classifier
Capper et al. (2018) developed a Random Forest classifier capable of differentiating between more than 80 CNS tumor classes and 9 control tissue samples based on methylation porfiles (data available in GEO repositories [GSE90496](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE90496) (training data) and [GSE109379](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109379) (test data)). This classifier was used as a basis for the current analysis. The pipeline and code used are described in further detail in the repository associated with the publication (Capper, 2018). This is provided by the authors upon request. The parts of their code used and adapted for the current analysis are labeled at the top of each script.

### Additional Data

Data from Brain Metastasis samples present in [GSE108576](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108576) and [GSE44661](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE44661) were used. Additionally, the methylation profiles of cell line LN229 were produced with Illumina EPIC Array and Nanopore Sequencing with Nanopolish (version 0.13.2) and Megalodon (version 2.2.10 with Guppy version 4.0.11 (See `./Data`).

### Preprocessing

Data preprocessing of .idat files from the previously mentioned GEO repositories was done with:

```r
code/preprocessing_MyDir.R
```
Data preprocessing of signal intensity files from Brain Metastasis data was done with:
```r
code/preprocessing_BrainMetastasis.R

```

### Training, Cross-Validation & Predictions

Random Forest training using methylation profiles was done with:

```r
code/training_MyDir.R
```
Random Forest training using the probes represented in a specific sample was done with:
```r
code/training_per_sample.R
```


Nested cross-validation was done with:

```r
code/cross-validation_MyDir.R
code/R/calculateCVfold.R
```
This performs the cross-validation on the modified Capper models and also performs the variable down-sampling experiments (See `./code/R/calculateCVfold.R`). Modify the file names in the script accordingly.

Training of calibration model was done with:

```r
code/calibration_MyDir.R
```

Prediction of methylation class of new samples was done with:

```r
code/R/predictionNewSamples.R
```

### tSNE analysis

tSNE visualizations were done with:

```r
code/tsne_mapping_MyDir.R
```

### Missing Value Imputation
Missing value imputation on Nanopore methylation calls followed three approaches: naive (replace `NA` with `0.5`), mean (replace `NA` with mean of probes from that sample), k-Approximate Nearest Neighbors (k-ANN) mean (replace `NA` with mean of k-ANN probes from that sample). This was done with:

```r
code/missing_value_imputation.R
code/R/na_impute.R
```

### Analysis of probes recovered by Nanopore
Identification of probes with missing values in Nanopore data, probes in the Illumina HM450 Array missing from Nanopore data, and probes in specified RF models missing from Nanopore data was done with:

```r
code/R/probe_analysis_newNanopore.R
```

### Performance Metrics

#### Misclassification error
Misclassification error of the predictions from cross-validation and predictions on training and test sets were estimated with:

```r
code/R/misclass_err_Capper_sets.R
code/R/misclass_err_CV.R
```

#### F1-score
Weighted F1-score (See `./code/R/f1_score_fold.R`) of the predictions from cross-validation and predictions on training and test sets were estimated with:

```r
code/R/f1_scores_Capper_sets.R
code/R/f1_scores_CV.R
```

#### Area Under Curve Precision-Recall (AUC-PR)
AUC-PR of the predictions from cross-validation were estimated with:
```r
code/PR-AUC_CV.R
```

### Additional scripts
Additional scripts are included. Among these `code/R/type2upper.R` allows to transform between methylation classes and methylation families, `code/R/waterfall_plots.R` produces waterfall plots from raw and calibrated class probabilities, and `code/R/confusionMatrices.R` generates and plots confusion matrices of the RF predictions.

## References
[Capper, D., Jones, D. T., Sill, M., Hovestadt, V., Schrimpf, D., Sturm, D., ... & Pfister, S. M. (2018). DNA methylation-based classification of central nervous system tumours. Nature, 555(7697), 469-474.](https://www.nature.com/articles/nature26000)
