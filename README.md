# Co-DEG & Machine Learning-Based Multi-Omics Analysis

This repository contains R scripts for integrated analysis of bulk microarray and single-cell RNA-seq data, focused on co-differentially expressed genes (Co-DEGs). The workflow includes functional interpretation (GO/KEGG), network analysis, GSVA, TF activity estimation, and predictive modeling using machine learning techniques.

## 🔧 Input Data

* **scRNA-seq Dataset**: GSE159677
* **Microarray Dataset**: GSE100927
* **Bulk Microarray DEG**: `DEG1_Disease_High_vs_Low.csv`, etc. (4 files total)
* **scRNA-seq DEG**: `scRNA_DEG1_*.csv`, etc. (4 files total)
* **Gene Expression Matrix**: `GSE100927_gene_mapped_expression.csv`

## 🧬 Co-DEG Comparison Definitions

Co-DEGs were identified by intersecting differentially expressed genes (DEGs) from microarray and single-cell RNA-seq analyses under the following conditions:

* **CoDEG1**: OASL High vs Low in **Disease** samples
* **CoDEG2**: OASL High vs Low in **NonDisease** samples
* **CoDEG3**: **Disease vs NonDisease** among **OASL High** samples
* **CoDEG4**: **Disease vs NonDisease** among **OASL Low** samples

## 📁 Workflow Overview

### 1. Co-DEG Identification

* Extract genes consistently differentially expressed in both bulk and scRNA-seq data
* Apply thresholds: |FC| > 1 and adjusted p < 0.05
* Annotate direction: UpUp, DownDown, Discordant

### 2. Functional Enrichment

* Perform GO and KEGG analysis using `clusterProfiler`
* Use only UpUp and DownDown genes
* Visualizations include barplots and radar charts

### 3. GSVA (Gene Set Variation Analysis)

* Based on MSigDB hallmark gene sets
* Split samples into OASL High vs Low expression groups
* Generate heatmaps, boxplots, and bar charts of pathway activities

### 4. Transcription Factor Activity

* Estimate TF activity using VIPER and dorothea regulons
* Visualize results via heatmaps and barplots

### 5. Machine Learning Models

* Apply LASSO, Random Forest, and SVM-RFE
* Select genes based on Co-DEG results
* Identify overlapping genes across models and evaluate prediction performance via ROC curves

### 6. ROC Analysis

* Evaluate predictive power of selected genes and model consensus genes
* Visualize ROC curves and AUC-based rankings

### 7. PPI Network Analysis

* Extract key GO-BP terms and build correlation-based pseudo PPI networks
* Highlight hub genes based on degree centrality

---

## 💻 Environment

* R (>= 4.2.0)
* Key packages:

  * `dplyr`, `readr`, `ggplot2`, `clusterProfiler`, `org.Hs.eg.db`
  * `GSVA`, `msigdbr`, `ComplexHeatmap`, `pheatmap`
  * `viper`, `dorothea`, `glmnet`, `randomForest`, `caret`, `pROC`

## 📂 Output Files

* `CoDEG*_Annotated.csv`: Annotated Co-DEG results
* `GO_Enrichment_*.png`, `KEGG_Enrichment_*.png`: Functional plots
* `Radar_*.png`: Radar plots of GO/KEGG enrichment
* `Volcano_*.png`: Volcano plots
* `TF Activity Heatmap`, `ROC Curve`, `AUC Barplot`, `PPI Network`

## 🔁 Reproducibility Notes

* All CSV files should be in the working directory
* Gene symbols must be consistently uppercase or lowercase
* The expression matrix file should contain gene symbols in the first column

## 📚 References
* GSE100927: [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100927]
* GSE159677: [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse159677]
* MSigDB: [https://www.gsea-msigdb.org/](https://www.gsea-msigdb.org/)
* Dorothea: [https://saezlab.github.io/dorothea/](https://saezlab.github.io/dorothea/)
* clusterProfiler: [https://doi.org/10.1093/bioinformatics/btq064](https://doi.org/10.1093/bioinformatics/btq064)

---

This pipeline allows integrated biological interpretation and machine learning-based evaluation of gene expression signatures across bulk and single-cell transcriptomics.
