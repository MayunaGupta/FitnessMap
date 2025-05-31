# GLS Analysis for Gene Interactions

## Overview

This repository contains code for performing Generalized Least Squares (GLS) regression analysis on gene expression data. The analysis identifies significant gene interactions and calculates fitness metrics based on GLS coefficients.

## Contents

- `GLS.py`: Main script for performing GLS analysis.
- `visualise_GLS.ipynb`: Jupyter notebook for visualizing the results of the GLS analysis.

## Features

- **GLS Regression**: Computes GLS coefficients and p-values for gene interactions.
- **FDR Correction**: Applies Benjamini-Hochberg FDR correction to control for multiple testing.
- **Visualization**: Generates heatmaps, histograms, and scatter plots to visualize significant interactions and fitness metrics.

## Installation

Make sure you have the following Python packages installed:

```bash
pip install pandas numpy scipy statsmodels matplotlib seaborn
```

## Usage

1. **Prepare Your Data**: Ensure your gene expression data is in a tab-separated values (TSV) format. The first few columns should contain metadata, while the remaining columns should contain numeric values representing gene expressions.

2. **Run the GLS Analysis**:
   - Open a terminal and navigate to the directory containing `GLS.py`.
   - Run the script using Python:
     ```bash
     python GLS.py
     ```

3. **Output Files**:
   - `ecoli_p_values.csv`: Contains the original p-values for gene interactions.
   - `ecoli_significant_interactions.csv`: Contains significant gene pairs after FDR correction.
   - `ecoli_gls_coefficients.csv`: Contains the GLS coefficients for all gene pairs.
   - `ecoli_fitness_results.csv`: Contains fitness metrics for significant gene pairs.

4. **Visualize Results**:
   - Open `visualise_GLS.ipynb` in Jupyter Notebook.
   - Run the cells to visualize the distribution of p-values, significant interactions, and fitness metrics.

## Analyzing Results

- **P-Value Distribution**: The histogram shows the distribution of original p-values, with a significance threshold marked at 0.05.
- **Top Interactions**: The top 10 most similar gene interactions based on GLS coefficients can be visualized in a heatmap.
- **Fitness Metrics**: Fitness metrics for significant gene pairs are calculated and can be analyzed further.

## Conclusion

This analysis provides insights into gene interactions and their significance based on GLS regression. The visualizations help in understanding the relationships between genes and their effects.

For any questions or issues, please feel free to reach out.