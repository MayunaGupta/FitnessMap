import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from statsmodels.regression.linear_model import GLS

# Load the fitness values file
df = pd.read_csv("fit_organism_Keio.tsv", sep="\t")

# Display the first few rows of the dataframe
print(df.head())

# Filter the dataframe to include only the columns of interest
# Select the 'geneName' column and any columns that contain 'set'
gene_columns = ['geneName'] + [col for col in df.columns if 'set' in col]
filtered_df = df[gene_columns]

# Sort the dataframe alphabetically by gene names and take the first 500 rows
# filtered_df = filtered_df.sort_values(by='geneName').head(500)

# Pivot the dataframe to create a matrix of genes vs. fitness values
fitness_matrix = filtered_df.set_index('geneName')

# Calculate the GLS correlation coefficients between genes
def gls_correlation(x, y):
    # Add constant term for regression
    X = np.column_stack([np.ones(len(x)), x])
    # Fit GLS model
    model = GLS(y, X)
    results = model.fit()
    # Calculate correlation from R-squared
    return np.sqrt(results.rsquared)

# Transpose the matrix to get genes as rows
fitness_matrix_t = fitness_matrix

# Create correlation matrix using GLS between genes
correlation_matrix = pd.DataFrame(index=fitness_matrix_t.columns, columns=fitness_matrix_t.columns)
# for i in range(len(fitness_matrix_t.columns)):
#     for j in range(len(fitness_matrix_t.columns)):
#         if i != j:
#             correlation_matrix.iloc[i, j] = gls_correlation(fitness_matrix_t.iloc[:, i], fitness_matrix_t.iloc[:, j])
#         else:
#             correlation_matrix.iloc[i, j] = 1.0

correlation_matrix = fitness_matrix_t.corr(method='pearson')


# Display the correlation matrix
print("\nCorrelation Matrix Shape:", correlation_matrix.shape)
print("\nFirst few rows and columns of correlation matrix:")
print(correlation_matrix.iloc[:5, :5])

# Save correlation matrix to a file
# correlation_matrix.to_csv('gene_correlations_gls.tsv', sep='\t')

# Display summary statistics of correlations
print("\nCorrelation Statistics:")
print("Mean correlation:", correlation_matrix.values[np.triu_indices_from(correlation_matrix, k=1)].mean())
print("Max correlation:", correlation_matrix.values[np.triu_indices_from(correlation_matrix, k=1)].max())
print("Min correlation:", correlation_matrix.values[np.triu_indices_from(correlation_matrix, k=1)].min())

# Identify the most correlated gene pairs
# Set a threshold for correlation
threshold = 0.8
highly_correlated_pairs = correlation_matrix[(correlation_matrix > threshold) & (correlation_matrix < 1.0)]

print("\nHighly correlated gene pairs (correlation >", threshold, "):")
print(highly_correlated_pairs)

# Save highly correlated pairs to a file
highly_correlated_pairs.to_csv('highly_correlated_genes_gls.tsv', sep='\t')

# Visualize the correlation matrix using a heatmap
plt.figure(figsize=(12, 10))
sns.heatmap(correlation_matrix, cmap='coolwarm', annot=True, fmt=".2f", square=True)
plt.title('Gene Correlation Matrix')
plt.show()
'''
# Optional: Create a network graph of highly correlated genes
import networkx as nx

# Create a graph from the correlation matrix
G = nx.from_pandas_adjacency(highly_correlated_pairs)

# Draw the network
plt.figure(figsize=(12, 12))
pos = nx.spring_layout(G)
nx.draw(G, pos, with_labels=True, node_size=700, node_color='lightblue', font_size=10, font_color='black', font_weight='bold', edge_color='gray')
plt.title('Network of Highly Correlated Genes')
plt.show()


# are there distinct clusters of highly correlated genes?
from scipy.cluster.hierarchy import linkage, dendrogram

# Perform hierarchical clustering
linked = linkage(correlation_matrix, method='ward')

# Create a dendrogram to visualize clusters
plt.figure(figsize=(12, 8))
dendrogram(linked, labels=correlation_matrix.index, orientation='top', leaf_rotation=90)
plt.title('Hierarchical Clustering of Genes')
plt.xlabel('Genes')
plt.ylabel('Distance')
plt.show()


# Calculate the degree of each node in the graph
degree_dict = dict(G.degree())
degree_values = list(degree_dict.values())

# Plot the degree distribution
plt.figure(figsize=(10, 6))
plt.hist(degree_values, bins=range(1, max(degree_values) + 1), alpha=0.75)
plt.title('Degree Distribution of Gene Correlation Network')
plt.xlabel('Degree')
plt.ylabel('Frequency')
plt.show()


# Define a threshold for hub genes
hub_threshold = 5  # Adjust this threshold as needed
hub_genes = [gene for gene, degree in degree_dict.items() if degree > hub_threshold]

# Display the hub genes
print("Hub Genes:", hub_genes)

'''