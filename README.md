# FitnessMap
Working on fitness map project
## Exploratory Data Analysis Questions

### Subset: first 500 pairs alphabettically for the Ecoli fitness: fit_organism_Keio.tsv file 

### 1. Distribution Analysis (Shay)
- What is the overall distribution of correlation coefficients? 
- Are correlations normally distributed?
- What are the ranges of correlation values?
- What proportion of correlations are positive vs negative?

### 2. Network Properties (Mayuna)
- Which genes have the most correlations with other genes?
- Are there distinct clusters of highly correlated genes?
- What is the degree distribution of the gene correlation network?
- Are there hub genes that correlate strongly with many others?

### 3. Functional Analysis (Kate)
- Are genes with similar functions more likely to be correlated?
- How do correlations differ between:
  * Essential vs non-essential genes
  * Different functional categories (metabolism, regulation, etc.)
  * Conserved vs non-conserved correlations

### 4. Pathway Analysis
- Do genes in the same metabolic pathway show stronger correlations?
- Are there unexpected correlations between different pathways?
- How do regulatory genes correlate with their target genes?

### 5. Condition-Specific Patterns
- Looking at the fitness data (specific_phenotypes_Keio.txt):
  * Which conditions show the strongest gene-gene correlations?
  * Are there condition-specific correlation patterns?
  * How do stress responses affect gene correlations?

### 6. Statistical Properties (Tristan)
- What is the significance threshold for correlations?
- How many correlations are statistically significant?
- What is the false discovery rate?
- Are there any systematic biases in the correlation patterns?

### 7. Comparative Analysis
- How do conserved correlations differ from non-conserved ones?
- Are there any evolutionary patterns in the correlations?
- Do essential genes show different correlation patterns?

### 8. Structural Analysis
- Is there a relationship between gene proximity on the chromosome and correlation strength?
- Do operons show distinct correlation patterns?
- How do regulatory elements affect correlation patterns?

### 9. Biological Validation
- Do known genetic interactions correspond to strong correlations?
- Can we validate correlations against known biological pathways?
- Are there novel correlations that suggest unknown interactions?

### 10. Technical Quality Assessment (Sarah)
- Are there any systematic biases in the data?
- How complete is the correlation network?
- Are there genes with missing or unreliable data?

## Clustering Analysis by Mayuna

### Overview
Developed a comprehensive functional clustering analysis system that groups genes by biological function and analyzes fitness patterns within each category. The analysis focuses on understanding how genes respond to different stress conditions and nutrient sources.

### Functional Categories Identified
The analysis categorizes the 173 experimental conditions into:
- **Carbon Sources** (marked with "(C)"): Alternative carbon metabolism conditions
- **Nitrogen Sources** (marked with "(N)"): Alternative nitrogen metabolism conditions  
- **Antibiotics**: Various antibiotic stress conditions (chloramphenicol, tetracycline, etc.)
- **Metal Stress**: Heavy metal toxicity conditions (nickel, cobalt, copper, etc.)
- **Chemical Stress**: Organic compound stress (DMSO, methylglyoxal, furfuraldehyde, etc.)
- **Motility**: Bacterial motility assays
- **Other**: Miscellaneous conditions

### Analysis Methods
1. **Data Processing**: 
   - Loaded fitness data from `fit_organism_Keio.tsv` (3789 genes × 173 conditions)
   - Filtered and cleaned data, removing genes with missing names
   - Transposed data for gene-wise analysis

2. **Functional Categorization**:
   - Automatically categorized experimental conditions by biological function
   - Created separate datasets for each functional category

3. **Clustering Analysis**:
   - Applied K-means clustering within each functional category
   - Used standardized fitness values for clustering
   - Generated 5 clusters per category by default

4. **Similarity Analysis**:
   - Calculated cosine similarity between gene pairs within categories
   - Identified high-fitness gene pairs (similarity > 0.7)
   - Analyzed combined fitness scores for gene pairs

5. **Correlation Analysis**:
   - Computed Pearson correlations between top genes within categories
   - Identified highly correlated gene pairs (|correlation| > 0.8)
   - Generated correlation matrices for network analysis

### Key Findings
- **Top Performing Genes**: Identified genes with highest average fitness in each category
- **Gene Pairs**: Found gene pairs with high combined fitness and similarity
- **Functional Patterns**: Revealed category-specific fitness patterns
- **Gene Networks**: Discovered highly correlated gene clusters within functional categories

### Files Created
- `functional_clustering.py`: Main analysis framework with FunctionalGeneClusterAnalyzer class
- `targeted_analysis_simple.py`: Simplified analysis script for specific categories
- `explore_columns.py`: Column categorization and exploration script
- `results_*/`: Category-specific result directories containing:
  - `top_genes.csv`: Top genes by average fitness
  - `high_fitness_pairs.csv`: Gene pairs with high combined fitness
  - `cluster_assignments.csv`: Cluster assignments for all genes
  - `gene_correlations.csv`: Correlation matrices
  - `summary_stats.csv`: Category summary statistics

### How to Generate Visualizations

#### 1. Run Complete Functional Analysis
```bash
# Activate your conda environment
conda activate QML

# Run the main functional clustering analysis
python functional_clustering.py
```
This generates:
- Interactive UMAP plots for each functional category
- Cluster visualizations with color-coded gene groups
- Automatic display of plots in browser

#### 2. Run Targeted Analysis
```bash
# Run simplified analysis for specific categories
python targeted_analysis_simple.py
```
This creates:
- Detailed analysis for carbon sources, nitrogen sources, antibiotics, and metal stress
- Category-specific result directories
- Summary comparison across categories
- Gene correlation analysis

#### 3. Generate Custom Visualizations
```python
# Example: Analyze specific category
from functional_clustering import FunctionalGeneClusterAnalyzer, categorize_columns
import pandas as pd

# Load and process data
df = pd.read_csv("fit_organism_Keio.tsv", sep="\t")
# ... data processing steps ...

# Create analyzer
analyzer = FunctionalGeneClusterAnalyzer(fitness_data, categories)

# Generate UMAP visualization for carbon sources
fig = analyzer.visualize_category_clusters('carbon_sources')
fig.show()

# Analyze gene networks
results, analyzer = analyze_specific_categories(['carbon_sources'])
analyze_gene_correlations(results, 'carbon_sources', top_n=20)
```

#### 4. View Results
- **Interactive plots**: Automatically open in web browser during analysis
- **CSV files**: Located in `results_*/` directories for each category
- **Summary files**: 
  - `category_summary_comparison.csv`: Cross-category comparison
  - `functional_insights.txt`: Key findings summary

### Dependencies Required
```bash
conda install scikit-learn umap-learn plotly pandas numpy
conda install -c conda-forge networkx xgboost
```

### Analysis Outputs
- **Cluster Assignments**: Genes grouped by fitness patterns within categories
- **Top Gene Lists**: Highest-performing genes in each functional category
- **Gene Pair Analysis**: Synergistic gene combinations with high fitness
- **Correlation Networks**: Gene interaction patterns within categories
- **Interactive Visualizations**: UMAP embeddings with cluster coloring

## Data Files
- `cofit_organism_Keio.txt`: Contains pre-calculated Pearson correlations between gene pairs
- `specific_phenotypes_Keio.txt`: Contains fitness data across different conditions
- `fit_organism_Keio.tsv`: Main fitness dataset (3789 genes × 173 conditions)
- Additional data files to be documented...

## Analysis Steps
1. Load and clean the data
2. Calculate basic statistics
3. Create visualizations for distributions and patterns
4. Cross-reference with other available data sources
