import pandas as pd
import numpy as np
from scipy.special import stdtr
from statsmodels.stats.multitest import fdrcorrection
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

class GLSRegression:
    def __init__(self, screens):
        # Convert screens array to a DataFrame and fill NaN with 0
        self.screens_df = pd.DataFrame(screens).fillna(0)
        self.Warped_Screens = None
        self.Warped_Intercept = None
        self.GLS_coef = None
        self.GLS_se = None
        self.GLS_pW_df = None

    def whiten_data(self):
        # Compute Cholesky decomposition of inverse covariance matrix
        # Note: We need to transpose the data to get the correct covariance matrix
        cholsigmainv = np.linalg.cholesky(np.linalg.inv(np.cov(self.screens_df.T)))
        # Apply the transformation to "whiten" the data
        self.Warped_Screens = self.screens_df.values @ cholsigmainv
        self.Warped_Intercept = cholsigmainv.sum(axis=0)

    def linear_regression(self):
        num_genes = self.Warped_Screens.shape[0]
        self.GLS_coef = np.zeros((num_genes, num_genes))
        self.GLS_se = np.zeros((num_genes, num_genes))
        ys = self.Warped_Screens.T  # Each column is a gene

        for gene_index in range(num_genes):
            # Set up design matrix: intercept + predictor gene
            x = np.stack((np.ones(len(self.Warped_Screens[gene_index])), self.Warped_Screens[gene_index]), axis=1)
            coef, residues, rank, s = np.linalg.lstsq(x, ys, rcond=None)

            # Fallback for zero residuals
            if residues.size == 0:
                residues = np.array([0])

            df = len(x) - x.shape[1]
            self.GLS_coef[gene_index, :] = coef[1, :]  # Slope = effect size
            se_matrix = np.sqrt(np.linalg.pinv(x.T @ x)[1, 1] * residues / df)
            self.GLS_se[gene_index, :] = se_matrix

    def calculate_p_values(self):
        # Calculate t-statistics
        t_stat = self.GLS_coef / self.GLS_se
        # Degrees of freedom: number of screens - 2 (accounts for intercept and slope)
        df = self.screens_df.shape[1] - 2
        # Compute two-tailed p-values using the Student's t-distribution
        GLS_pW = 2 * (1 - stdtr(df, np.abs(t_stat)))
        # Set diagonal values to 1 (no self-comparison)
        np.fill_diagonal(GLS_pW, 1)
        
        # Create index for the DataFrame
        gene_indices = [f"Gene_{i}" for i in range(GLS_pW.shape[0])]
        
        # Store the result in a DataFrame for easier handling
        self.GLS_pW_df = pd.DataFrame(GLS_pW, 
                                     index=gene_indices, 
                                     columns=gene_indices)

    def apply_fdr_correction(self, significance_threshold=0.005):
        # Stack the upper triangle of the p-value matrix
        stacked_p = self.GLS_pW_df.stack()
        stacked_p = stacked_p[stacked_p.index.get_level_values(0) < stacked_p.index.get_level_values(1)]
        # Apply Benjamini-Hochberg FDR correction
        fdr_corrected = pd.Series(fdrcorrection(stacked_p)[1], index=stacked_p.index)
        # Filter for gene pairs below threshold
        significant_p_values = fdr_corrected[fdr_corrected < significance_threshold]
        return significant_p_values

def analyze_ecoli_data():
    # Read the E. coli data
    ecoli_data = pd.read_csv('fit_organism_Keio.tsv', sep='\t')
    
    # Print the shape and column names to debug
    print("Original data shape:", ecoli_data.shape)
    print("Column names:", ecoli_data.columns.tolist())
    
    # Get only the numeric columns, excluding any metadata columns
    numeric_cols = ecoli_data.select_dtypes(include=[np.number]).columns
    print("Number of numeric columns:", len(numeric_cols))
    
    # Create a new DataFrame with only numeric columns
    screens_data = ecoli_data[numeric_cols]
    
    # Print the shape of the processed data
    print("Processed data shape:", screens_data.shape)
    
    # Run GLS analysis
    gls_model = GLSRegression(screens_data)
    gls_model.whiten_data()
    gls_model.linear_regression()
    gls_model.calculate_p_values()
    gls_model.GLS_pW_df.to_csv('ecoli_p_values.csv')
    significant_results = gls_model.apply_fdr_correction()
    significant_results = significant_results.reset_index()
    significant_results.columns = ['gene1', 'gene2', 'fdr_corrected_p_value']
    significant_results.to_csv('ecoli_significant_interactions.csv')
    
    # Create visualizations
    # 1. Heatmap of significant interactions
    # fig1, ax1 = plt.subplots(figsize=(12, 10))
    # significant_matrix = pd.DataFrame(0, 
    #                                index=gls_model.GLS_pW_df.index,
    #                                columns=gls_model.GLS_pW_df.columns)
    # for idx in significant_results.index:
    #     significant_matrix.loc[idx[0], idx[1]] = 1
    #     significant_matrix.loc[idx[1], idx[0]] = 1
    
    # sns.heatmap(significant_matrix, cmap='YlOrRd', ax=ax1)
    # ax1.set_title('Significant Gene-Gene Interactions')
    # plt.savefig('ecoli_heatmap.pdf')
    # plt.close()
    
    # # 2. Distribution of effect sizes
    # fig2, ax2 = plt.subplots(figsize=(10, 6))
    # sns.histplot(gls_model.GLS_coef.flatten(), bins=50, ax=ax2)
    # ax2.set_title('Distribution of Effect Sizes')
    # ax2.set_xlabel('Effect Size')
    # ax2.set_ylabel('Frequency')
    # plt.savefig('ecoli_effect_sizes.pdf')
    # plt.close()
    
    # 3. Volcano plot
    # fig3, ax3 = plt.subplots(figsize=(10, 8))
    # ax3.scatter(gls_model.GLS_coef.flatten(), 
    #            -np.log10(gls_model.GLS_pW_df.values.flatten()),
    #            alpha=0.5)
    # ax3.axhline(y=-np.log10(0.05), color='r', linestyle='--')
    # ax3.set_title('Volcano Plot')
    # ax3.set_xlabel('Effect Size')
    # ax3.set_ylabel('-log10(p-value)')
    # plt.savefig('ecoli_volcano.pdf')
    # plt.close()
    
    return significant_results


def analyze_ecoli_data2():
    # Read the E. coli data
    ecoli_data = pd.read_csv('fit_organism_Keio.tsv', sep='\t')
    
    # Print the shape and column names to debug
    print("Original data shape:", ecoli_data.shape)
    print("Column names:", ecoli_data.columns.tolist())
    
    # Get only the numeric columns, excluding any metadata columns
    numeric_cols = ecoli_data.select_dtypes(include=[np.number]).columns
    print("Number of numeric columns:", len(numeric_cols))
    
    # Create a new DataFrame with only numeric columns
    screens_data = ecoli_data[numeric_cols]
    
    # Print the shape of the processed data
    print("Processed data shape:", screens_data.shape)
    
    # Run GLS analysis
    gls_model = GLSRegression(screens_data)
    gls_model.whiten_data()
    gls_model.linear_regression()
    gls_model.calculate_p_values()
    gls_model.GLS_pW_df.to_csv('ecoli_p_values.csv')
    significant_results = gls_model.apply_fdr_correction()
    significant_results = significant_results.reset_index()
    significant_results.columns = ['gene1', 'gene2', 'fdr_corrected_p_value']
    significant_results.to_csv('ecoli_significant_interactions.csv')

    # Analyze fitness of gene pairs based on GLS coefficients
    fitness_results = analyze_fitness(gls_model.GLS_coef, significant_results)
    fitness_results.to_csv('ecoli_fitness_results.csv')

    # Save the GLS coefficients to a CSV file
    coef_df = pd.DataFrame(gls_model.GLS_coef, index=[f"Gene_{i}" for i in range(gls_model.GLS_coef.shape[0])])
    coef_df.to_csv('ecoli_gls_coefficients.csv')

    return significant_results

def analyze_fitness(GLS_coef, significant_results):
    # Create a DataFrame to hold fitness results
    fitness_data = []

    for index, row in significant_results.iterrows():
        gene1_index = int(row['gene1'].split('_')[1])  # Extract gene index from name
        gene2_index = int(row['gene2'].split('_')[1])  # Extract gene index from name
        
        # Calculate fitness metric (e.g., mean of coefficients for the gene pair)
        fitness_metric = np.mean([GLS_coef[gene1_index, gene2_index], GLS_coef[gene2_index, gene1_index]])
        
        fitness_data.append({
            'gene1': row['gene1'],
            'gene2': row['gene2'],
            'fitness_metric': fitness_metric
        })

    return pd.DataFrame(fitness_data)

if __name__ == "__main__":
    results = analyze_ecoli_data2()
    print(f"Number of significant interactions: {len(results)}")
# if __name__ == "__main__":
#     results = analyze_ecoli_data()
#     print(f"Number of significant interactions: {len(results)}")

#     results.to_csv('ecoli_significant_interactions.csv')

