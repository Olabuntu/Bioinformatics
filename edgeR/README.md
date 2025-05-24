This project performs differential gene expression analysis for RNA-seq and Ribo-seq data using edgeR. The analysis includes data preprocessing, normalization, quality control, statistical modeling, and visualization of results.

## Steps in the Analysis

1. **Load Required Libraries**  
  The script loads necessary R libraries for data manipulation, visualization, and statistical analysis.

2. **Load Data**  
  RNA-seq and Ribo-seq count data are loaded from CSV files.

3. **Data Preprocessing**  
  - Experimental groups are defined.
  - Lowly expressed genes are filtered out.
  - Data is normalized using the TMM method.

4. **Quality Control**  
  - MDS plots are generated to visualize sample relationships.
  - Log2 CPM boxplots are created to assess data distribution.

5. **Model Fitting and Dispersion Estimation**  
  - Design matrices are created for RNA-seq and Ribo-seq data.
  - Dispersion is estimated, and mean-dispersion relationships are visualized.

6. **Differential Expression Analysis**  
  - Generalized linear models are fitted.
  - Differential expression is tested, and results are saved to CSV files.

7. **Significant Results**  
  - Genes with FDR < 0.05 are identified and saved.
  - Volcano and MA plots are generated to visualize significant genes.

8. **Heatmap Visualization**  
  - Heatmaps are created for all genes and significant genes.

## Output Files
- Quality control plots: `MDS_rna.png`, `MDS_ribo.png`, `logCPM_rna.png`, `logCPM_ribo.png`
- Dispersion plots: `mean_dispersion_rna.png`, `mean_dispersion_ribo.png`
- Differential expression results: `rna_results.csv`, `ribo_results.csv`, `rna_results_sig.csv`, `ribo_results_sig.csv`
- Volcano plots: `volcano_rna.png`, `volcano_ribo.png`
- MA plots: `MA_rna_legend.png`, `MA_ribo_legend.png`
- Heatmaps: `heatmap_rna.png`, `heatmap_ribo.png`, `heatmap_rna_sig.png`, `heatmap_ribo_sig.png`

## Prerequisites
- Install required R libraries using `install.packages()` or `BiocManager::install()` for Bioconductor packages.
- Ensure input files `rna_counts.csv` and `ribo_counts.csv` are in the working directory.

## Usage
Run the script in R or RStudio. The results will be saved in the `new_results` directory.
