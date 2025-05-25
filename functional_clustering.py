import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.metrics.pairwise import cosine_similarity
import umap
import plotly.express as px
import plotly.graph_objects as go
from typing import Tuple, List, Dict
import xgboost as xgb
from sklearn.cluster import KMeans
import networkx as nx

def categorize_columns(columns):
    """
    Categorize columns by their functional groups.
    """
    categories = {
        'carbon_sources': [],
        'nitrogen_sources': [],
        'antibiotics': [],
        'metal_stress': [],
        'chemical_stress': [],
        'motility': [],
        'other': []
    }
    
    for col in columns:
        if 'set' not in col:
            continue
            
        col_lower = col.lower()
        
        # Carbon sources (marked with (C))
        if '(c)' in col_lower:
            categories['carbon_sources'].append(col)
        # Nitrogen sources (marked with (N))
        elif '(n)' in col_lower:
            categories['nitrogen_sources'].append(col)
        # Antibiotics
        elif any(antibiotic in col_lower for antibiotic in [
            'chloramphenicol', 'tetracycline', 'spectinomycin', 'carbenicillin',
            'cephalothin', 'bacitracin', 'fusidic', 'doxycycline', 'nalidixic',
            'phosphomycin', 'cycloserine'
        ]):
            categories['antibiotics'].append(col)
        # Metal stress
        elif any(metal in col_lower for metal in [
            'nickel', 'cobalt', 'copper', 'aluminum', 'thallium', 'cisplatin'
        ]):
            categories['metal_stress'].append(col)
        # Chemical stress
        elif any(chemical in col_lower for chemical in [
            'dmso', 'dimethyl sulfoxide', 'methylglyoxal', 'furfuraldehyde',
            'hydroxymethylfurfural', 'vanillin', 'syringaldehyde', 'benzoic',
            'fluoride', 'chlorite', 'chloride'
        ]):
            categories['chemical_stress'].append(col)
        # Motility assays
        elif 'motility' in col_lower:
            categories['motility'].append(col)
        else:
            categories['other'].append(col)
    
    return categories

class FunctionalGeneClusterAnalyzer:
    def __init__(self, data: pd.DataFrame, functional_categories: Dict[str, List[str]]):
        """
        Initialize the analyzer with functional categories.
        
        Args:
            data: DataFrame with genes as rows and conditions as columns
            functional_categories: Dictionary mapping category names to column lists
        """
        self.data = data
        self.categories = functional_categories
        self.category_data = {}
        self.category_clusters = {}
        
        # Prepare data for each category
        for category, columns in functional_categories.items():
            if columns:  # Only process categories with columns
                self.category_data[category] = data[columns]
    
    def cluster_by_category(self, category: str, n_clusters: int = 5) -> np.ndarray:
        """
        Perform clustering for a specific functional category.
        
        Args:
            category: Name of the functional category
            n_clusters: Number of clusters
            
        Returns:
            Array of cluster labels
        """
        if category not in self.category_data:
            print(f"Category '{category}' not found or has no data")
            return np.array([])
        
        data = self.category_data[category]
        scaler = StandardScaler()
        scaled_data = scaler.fit_transform(data)
        
        # Use KMeans for clustering
        kmeans = KMeans(n_clusters=n_clusters, random_state=42)
        clusters = kmeans.fit_predict(scaled_data)
        
        self.category_clusters[category] = clusters
        return clusters
    
    def visualize_category_clusters(self, category: str) -> go.Figure:
        """
        Create UMAP visualization for a specific category.
        
        Args:
            category: Name of the functional category
            
        Returns:
            Interactive UMAP plot
        """
        if category not in self.category_data or category not in self.category_clusters:
            print(f"No data or clusters for category '{category}'")
            return go.Figure()
        
        data = self.category_data[category]
        clusters = self.category_clusters[category]
        
        # Scale data
        scaler = StandardScaler()
        scaled_data = scaler.fit_transform(data)
        
        # Create UMAP embedding
        reducer = umap.UMAP(
            n_neighbors=min(15, len(data)-1),
            min_dist=0.1,
            n_components=2,
            random_state=42
        )
        embedding = reducer.fit_transform(scaled_data)
        
        # Create plot
        df_plot = pd.DataFrame({
            'UMAP1': embedding[:, 0],
            'UMAP2': embedding[:, 1],
            'Cluster': clusters,
            'Gene': list(data.index)
        })
        
        fig = px.scatter(
            df_plot,
            x='UMAP1',
            y='UMAP2',
            color='Cluster',
            hover_data=['Gene'],
            title=f'UMAP Visualization - {category.replace("_", " ").title()}',
            template='plotly_white'
        )
        
        return fig
    
    def analyze_fitness_by_category(self, category: str, threshold: float = 0.8) -> Dict:
        """
        Analyze gene pairs within a functional category for high fitness.
        
        Args:
            category: Name of the functional category
            threshold: Minimum cosine similarity threshold
            
        Returns:
            Dictionary with fitness analysis results
        """
        if category not in self.category_data:
            return {}
        
        data = self.category_data[category]
        scaler = StandardScaler()
        scaled_data = scaler.fit_transform(data)
        
        # Calculate cosine similarity
        cosine_sim = cosine_similarity(scaled_data)
        
        results = {
            'high_fitness_pairs': [],
            'category': category,
            'n_conditions': len(data.columns),
            'mean_fitness': data.mean().mean(),
            'std_fitness': data.mean().std()
        }
        
        # Find high-fitness gene pairs
        for i, gene1 in enumerate(data.index):
            for j, gene2 in enumerate(data.index[i+1:], i+1):
                similarity = cosine_sim[i, j]
                
                if similarity >= threshold:
                    gene1_fitness = data.iloc[i]
                    gene2_fitness = data.iloc[j]
                    combined_fitness = (gene1_fitness + gene2_fitness) / 2
                    
                    if combined_fitness.mean() > results['mean_fitness']:
                        results['high_fitness_pairs'].append({
                            'gene1': gene1,
                            'gene2': gene2,
                            'similarity': similarity,
                            'combined_fitness': combined_fitness.mean(),
                            'fitness_std': combined_fitness.std()
                        })
        
        # Sort by combined fitness
        results['high_fitness_pairs'].sort(key=lambda x: x['combined_fitness'], reverse=True)
        
        return results
    
    def get_top_genes_by_category(self, category: str, n: int = 20) -> pd.DataFrame:
        """
        Get top genes by average fitness in a category.
        
        Args:
            category: Name of the functional category
            n: Number of top genes to return
            
        Returns:
            DataFrame with top genes and their fitness values
        """
        if category not in self.category_data:
            return pd.DataFrame()
        
        data = self.category_data[category]
        gene_fitness = data.mean(axis=1).sort_values(ascending=False)
        
        top_genes = gene_fitness.head(n)
        
        return pd.DataFrame({
            'Gene': top_genes.index,
            'Average_Fitness': top_genes.values,
            'Category': category
        })

def run_functional_analysis(data_file: str = "fit_organism_Keio.tsv"):
    """
    Run the complete functional clustering analysis.
    
    Args:
        data_file: Path to the fitness data file
    """
    # Load data
    df = pd.read_csv(data_file, sep="\t")
    
    # Get gene columns and fitness columns
    gene_columns = ['geneName'] + [col for col in df.columns if 'set' in col]
    filtered_df = df[gene_columns]
    
    # Remove rows with missing gene names
    filtered_df = filtered_df.dropna(subset=['geneName'])
    
    # Set gene names as index
    fitness_data = filtered_df.set_index('geneName')
    
    # Categorize columns
    fitness_columns = [col for col in fitness_data.columns if 'set' in col]
    categories = categorize_columns(fitness_columns)
    
    print("Functional Categories Found:")
    for category, columns in categories.items():
        print(f"  {category}: {len(columns)} conditions")
    
    # Initialize analyzer
    analyzer = FunctionalGeneClusterAnalyzer(fitness_data, categories)
    
    # Analyze each category
    results = {}
    for category in categories.keys():
        if category in analyzer.category_data and len(analyzer.category_data[category].columns) > 0:
            print(f"\nAnalyzing {category}...")
            
            # Cluster genes
            clusters = analyzer.cluster_by_category(category, n_clusters=5)
            
            # Create visualization
            fig = analyzer.visualize_category_clusters(category)
            fig.show()
            
            # Analyze fitness
            fitness_results = analyzer.analyze_fitness_by_category(category, threshold=0.7)
            
            # Get top genes
            top_genes = analyzer.get_top_genes_by_category(category, n=10)
            
            results[category] = {
                'clusters': clusters,
                'fitness_analysis': fitness_results,
                'top_genes': top_genes
            }
            
            print(f"Top 5 high-fitness gene pairs in {category}:")
            for i, pair in enumerate(fitness_results['high_fitness_pairs'][:5]):
                print(f"  {i+1}. {pair['gene1']} - {pair['gene2']}: {pair['combined_fitness']:.3f}")
    
    return results, analyzer

if __name__ == "__main__":
    results, analyzer = run_functional_analysis() 