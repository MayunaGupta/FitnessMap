# Gene Interaction Analysis Pipeline

## Overview

This repository contains a comprehensive pipeline for analyzing gene interactions in E. coli using Generalized Least Squares (GLS) regression and AlphaFold2 protein structure prediction. The analysis identifies significant gene interactions from fitness data and validates them through protein sequence and structure analysis.

## Project Structure

```
FitnessMap/
├── GLS.py                                    # Main GLS analysis script
├── visualise_GLS.ipynb                       # Visualization notebook
├── 250601_Get_protein_sequences_for_gene_pairs.ipynb    # Protein sequence extraction
├── 250601_AlphaFold2_Protein_Interactions.ipynb         # AlphaFold2 predictions
├── 250601_Get_min_PAE_interaction_Alphafold_prediction.ipynb  # PAE analysis
└── README.md                                 # This file
```

## Pipeline Overview

### 1. GLS Analysis (`GLS.py`)
- **Purpose**: Identifies significant gene interactions using GLS regression
- **Input**: `fit_organism_Keio.tsv` (fitness data)
- **Output**: 
  - `ecoli_p_values.csv` (p-values for all gene pairs)
  - `ecoli_significant_interactions.csv` (significant interactions after FDR)
  - `ecoli_gls_coefficients.csv` (GLS coefficients)
  - `ecoli_fitness_results.csv` (fitness metrics)

### 2. Protein Sequence Analysis (`250601_Get_protein_sequences_for_gene_pairs.ipynb`)
- **Purpose**: Extracts protein sequences for gene pairs
- **Inputs**:
  - `ecoli_gls_coefficients.csv`
  - `fit_organism_Keio.tsv`
  - `organism_Keio.faa`
  - `uniprotkb_ecoli_AND_model_organism_8333_2025_06_02.tsv`
- **Output**: `top_interaction_sequences.csv`

### 3. AlphaFold2 Predictions (`250601_AlphaFold2_Protein_Interactions.ipynb`)
- **Purpose**: Predicts protein-protein interactions using AlphaFold2
- **Input**: `top_interaction_sequences.csv`
- **Output**: ZIP files containing AlphaFold2 predictions
- **Note**: Designed to run on Google Colab with GPU support

### 4. PAE Analysis (`250601_Get_min_PAE_interaction_Alphafold_prediction.ipynb`)
- **Purpose**: Analyzes prediction accuracy using PAE metrics
- **Inputs**:
  - `top_interaction_sequences.csv`
  - AlphaFold2 output ZIP files
- **Output**: Updated dataframe with PAE metrics
- **Interpretation**:
  - Min PAE < 1.7: Strong evidence for interaction
  - Min PAE < 10: Possible interaction

## Installation

```bash
pip install pandas numpy scipy statsmodels matplotlib seaborn
```

## Data Requirements

All required data files are available in the shared Google Drive folder:
https://drive.google.com/drive/folders/1-CfFsDSErSE719UJ20IgW3lo0sXQXNh1

Required files:
1. `fit_organism_Keio.tsv` - Fitness values from Fitness Browser
2. `organism_Keio.faa` - Protein sequences from Fitness Browser
3. `uniprotkb_ecoli_AND_model_organism_8333_2025_06_02.tsv` - UniProt E. coli data

## Usage

### 1. Run GLS Analysis
```bash
python GLS.py
```

### 2. Extract Protein Sequences
- Open `250601_Get_protein_sequences_for_gene_pairs.ipynb`
- Run all cells to generate sequence data

### 3. Run AlphaFold2 Predictions
- Open `250601_AlphaFold2_Protein_Interactions.ipynb` in Google Colab
- Modify the gene pair indices as needed
- Run predictions
- Upload results to shared folder: https://drive.google.com/drive/folders/1Gz7YlvsNBLRM48aVSyd6X7ODeR3TPklX

### 4. Analyze PAE Metrics
- Open `250601_Get_min_PAE_interaction_Alphafold_prediction.ipynb`
- Run analysis on AlphaFold2 outputs

## Key Features

### GLS Analysis
- Handles missing values and data normalization
- Accounts for experimental dependencies
- Applies FDR correction for multiple testing
- Generates comprehensive visualizations

### Protein Analysis
- Combines data from multiple sources
- Handles missing sequences (e.g., RNA-encoding genes)
- Integrates with AlphaFold2 for structure prediction
- Provides quantitative metrics for interaction confidence

## Results Interpretation

### Statistical Significance
- FDR-corrected p-values < 0.005 indicate significant interactions
- GLS coefficients range from -5.03 to +5.03
- Positive coefficients suggest cooperative interactions
- Negative coefficients suggest antagonistic interactions

### Structural Validation
- Min PAE < 1.7: High confidence in predicted interaction
- Min PAE < 10: Possible interaction, requires further validation
- Higher PAE values: Less confident predictions

## Performance Notes

- AlphaFold2 predictions are computationally intensive
- Current implementation uses T4 GPUs on Google Colab
- Multiple Colab instances can be run in parallel
- Consider using A100 GPUs for faster processing

## Future Improvements

1. Optimize AlphaFold2 prediction speed
2. Implement parallel processing for GLS analysis
3. Add more visualization options
4. Integrate with additional protein databases
5. Develop automated validation pipeline

## Contributing

Feel free to submit issues and enhancement requests!

## Contributors
1. Mayuna Gupta, Sarah Johnson, Tristan Ferdinand, Kate Zhou, Akshay Uppal, Adrian Jinich

## Contact

For questions or issues, please reach out to the repository maintainers.
