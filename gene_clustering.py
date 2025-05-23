import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.metrics.pairwise import cosine_similarity
import umap.umap_ as umap
import plotly.express as px
import plotly.graph_objects as go
from typing import Tuple, List, Dict
import xgboost as xgb
from sklearn.cluster import KMeans
import networkx as nx

class GeneClusterAnalyzer:
    def __init__(self, data: pd.DataFrame, n_clusters: int = 5):
        """
        Initialize the GeneClusterAnalyzer with data and number of clusters.
        
        Args:
            data: DataFrame with genes as rows and conditions as columns
            n_clusters: Number of clusters to create
        """
        self.data = data
        self.n_clusters = n_clusters
        self.scaler = StandardScaler()
        self.scaled_data = self.scaler.fit_transform(data)
        self.clusters = None
        self.umap_embedding = None
        
    def perform_xgboost_clustering(self) -> np.ndarray:
        """
        Perform clustering using XGBoost's built-in capabilities.
        
        Returns:
            Array of cluster labels
        """
        # Convert data to DMatrix format
        dtrain = xgb.DMatrix(self.scaled_data)
        
        # Set up parameters for clustering
        params = {
            'objective': 'multi:softmax',
            'num_class': self.n_clusters,
            'max_depth': 5,
            'eta': 0.1,
            'silent': 1
        }
        
        # Train the model
        num_round = 100
        bst = xgb.train(params, dtrain, num_round)
        
        # Get predictions (cluster assignments)
        preds = bst.predict(dtrain)
        self.clusters = preds.astype(int)
        
        return self.clusters
    
    def perform_random_forest_clustering(self) -> np.ndarray:
        """
        Perform clustering using Random Forest and KMeans.
        
        Returns:
            Array of cluster labels
        """
        # Train Random Forest
        rf = RandomForestClassifier(
            n_estimators=100,
            max_depth=5,
            random_state=42
        )
        
        # Use Random Forest features for clustering
        rf_features = rf.fit_transform(self.scaled_data)
        
        # Apply KMeans clustering
        kmeans = KMeans(n_clusters=self.n_clusters, random_state=42)
        self.clusters = kmeans.fit_predict(rf_features)
        return self.clusters
    
    def create_umap_visualization(self) -> Tuple[np.ndarray, go.Figure]:
        """
        Create UMAP visualization of the data.
        
        Returns:
            Tuple of (UMAP embedding, interactive plot)
        """
        # Create UMAP embedding
        reducer = umap.UMAP(
            n_neighbors=15,
            min_dist=0.1,
            n_components=2,
            random_state=42
        )
        self.umap_embedding = reducer.fit_transform(self.scaled_data)
        
        # Create interactive plot
        df_plot = pd.DataFrame({
            'UMAP1': self.umap_embedding[:, 0],
            'UMAP2': self.umap_embedding[:, 1],
            'Cluster': self.clusters,
            'Gene': self.data.index
        })
        
        fig = px.scatter(
            df_plot,
            x='UMAP1',
            y='UMAP2',
            color='Cluster',
            hover_data=['Gene'],
            title='UMAP Visualization of Gene Clusters'
        )
        
        return self.umap_embedding, fig
    
    def analyze_cosine_distances(self, threshold: float = 0.8) -> Dict[int, List[Tuple[str, str, float]]]:
        """
        Analyze cosine distances between genes within each cluster.
        
        Args:
            threshold: Minimum cosine similarity threshold
            
        Returns:
            Dictionary mapping cluster numbers to lists of gene pairs and their similarities
        """
        cosine_sim = cosine_similarity(self.scaled_data)
        cluster_pairs = {}
        
        for cluster in range(self.n_clusters):
            cluster_genes = self.data.index[self.clusters == cluster]
            cluster_pairs[cluster] = []
            
            for i, gene1 in enumerate(cluster_genes):
                for gene2 in cluster_genes[i+1:]:
                    idx1 = self.data.index.get_loc(gene1)
                    idx2 = self.data.index.get_loc(gene2)
                    similarity = cosine_sim[idx1, idx2]
                    
                    if similarity >= threshold:
                        cluster_pairs[cluster].append((gene1, gene2, similarity))
        
        return cluster_pairs
    
    def create_network_visualization(self, cluster_pairs: Dict[int, List[Tuple[str, str, float]]]) -> go.Figure:
        """
        Create interactive network visualization of gene pairs.
        
        Args:
            cluster_pairs: Dictionary of gene pairs and their similarities
            
        Returns:
            Interactive network plot
        """
        G = nx.Graph()
        
        # Add edges with weights
        for cluster, pairs in cluster_pairs.items():
            for gene1, gene2, similarity in pairs:
                G.add_edge(gene1, gene2, weight=similarity)
        
        # Create layout
        pos = nx.spring_layout(G)
        
        # Create edge trace
        edge_x = []
        edge_y = []
        edge_weights = []
        for edge in G.edges(data=True):
            x0, y0 = pos[edge[0]]
            x1, y1 = pos[edge[1]]
            edge_x.extend([x0, x1, None])
            edge_y.extend([y0, y1, None])
            edge_weights.append(edge[2]['weight'])
        
        edge_trace = go.Scatter(
            x=edge_x, y=edge_y,
            line=dict(width=0.5, color='#888'),
            hoverinfo='none',
            mode='lines'
        )
        
        # Create node trace
        node_x = []
        node_y = []
        node_text = []
        for node in G.nodes():
            x, y = pos[node]
            node_x.append(x)
            node_y.append(y)
            node_text.append(node)
        
        node_trace = go.Scatter(
            x=node_x, y=node_y,
            mode='markers+text',
            hoverinfo='text',
            text=node_text,
            textposition="top center",
            marker=dict(
                showscale=True,
                colorscale='YlGnBu',
                size=10,
                colorbar=dict(
                    thickness=15,
                    title='Node Connections',
                    xanchor='left',
                    titleside='right'
                )
            )
        )
        
        fig = go.Figure(data=[edge_trace, node_trace],
                       layout=go.Layout(
                           title='Gene Interaction Network',
                           showlegend=False,
                           hovermode='closest',
                           margin=dict(b=20,l=5,r=5,t=40)
                       ))
        
        return fig

def run_analysis(data: pd.DataFrame, n_clusters: int = 5, similarity_threshold: float = 0.8):
    """
    Run the complete gene clustering analysis.
    
    Args:
        data: DataFrame with genes as rows and conditions as columns
        n_clusters: Number of clusters to create
        similarity_threshold: Minimum cosine similarity threshold
    """
    # Initialize analyzer
    analyzer = GeneClusterAnalyzer(data, n_clusters)
    
    # Perform clustering (you can choose either method)
    print("Performing XGBoost clustering...")
    clusters = analyzer.perform_xgboost_clustering()
    
    # Create UMAP visualization
    print("Creating UMAP visualization...")
    _, umap_fig = analyzer.create_umap_visualization()
    umap_fig.show()
    
    # Analyze cosine distances
    print("Analyzing cosine distances...")
    cluster_pairs = analyzer.analyze_cosine_distances(similarity_threshold)
    
    # Print cluster pair statistics
    for cluster, pairs in cluster_pairs.items():
        print(f"\nCluster {cluster} has {len(pairs)} gene pairs with similarity > {similarity_threshold}")
        if pairs:
            print("Top 5 pairs:")
            for gene1, gene2, sim in sorted(pairs, key=lambda x: x[2], reverse=True)[:5]:
                print(f"{gene1} - {gene2}: {sim:.3f}")
    
    # Create network visualization
    print("\nCreating network visualization...")
    network_fig = analyzer.create_network_visualization(cluster_pairs)
    network_fig.show()

# Example usage:
if __name__ == "__main__":
    # Load your data
    data = pd.read_csv("gene_correlations_gls.tsv", sep='\t', index_col=0)
    
    # Run analysis
    run_analysis(data, n_clusters=5, similarity_threshold=0.8) 