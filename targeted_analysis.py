import pandas as pd
import numpy as np
from functional_clustering import FunctionalGeneClusterAnalyzer, categorize_columns
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.metrics.pairwise import cosine_similarity
import plotly.express as px
import plotly.graph_objects as go

def analyze_specific_categories(categories_to_analyze=['carbon_sources', 'nitrogen_sources', 'antibiotics']):
    """
    Perform detailed analysis on specific functional categories.
    
    Args:
        categories_to_analyze: List of category names to focus on
    """
    # Load data
    df = pd.read_csv("fit_organism_Keio.tsv", sep="\t")
    
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
    
    # Initialize analyzer
    analyzer = FunctionalGeneClusterAnalyzer(fitness_data, categories)
    
    results = {}
    
    for category in categories_to_analyze:
        if category in analyzer.category_data and len(analyzer.category_data[category].columns) > 0:
            print(f"\n{'='*60}")
            print(f"DETAILED ANALYSIS: {category.upper().replace('_', ' ')}")
            print(f"{'='*60}")
            
            data = analyzer.category_data[category]
            print(f"Number of conditions: {len(data.columns)}")
            print(f"Number of genes: {len(data.index)}")
            
            # Basic statistics
            print(f"\nFitness Statistics:")
            print(f"  Mean fitness: {data.mean().mean():.3f}")
            print(f"  Std fitness: {data.mean().std():.3f}")
            print(f"  Min fitness: {data.min().min():.3f}")
            print(f"  Max fitness: {data.max().max():.3f}")
            
            # Cluster genes
            clusters = analyzer.cluster_by_category(category, n_clusters=5)
            
            # Get top genes by average fitness
            top_genes = analyzer.get_top_genes_by_category(category, n=20)
            print(f"\nTop 10 genes by average fitness:")
            for i, row in top_genes.head(10).iterrows():
                print(f"  {i+1:2d}. {row['Gene']}: {row['Average_Fitness']:.3f}")
            
            # Analyze fitness pairs
            fitness_results = analyzer.analyze_fitness_by_category(category, threshold=0.7)
            
            print(f"\nTop 10 high-fitness gene pairs (similarity > 0.7):")
            for i, pair in enumerate(fitness_results['high_fitness_pairs'][:10]):
                print(f"  {i+1:2d}. {pair['gene1']} - {pair['gene2']}: "
                      f"fitness={pair['combined_fitness']:.3f}, "
                      f"similarity={pair['similarity']:.3f}")
            
            # Save detailed results
            save_category_results(category, data, top_genes, fitness_results, clusters)
            
            results[category] = {
                'data': data,
                'top_genes': top_genes,
                'fitness_results': fitness_results,
                'clusters': clusters
            }
    
    return results, analyzer

def save_category_results(category, data, top_genes, fitness_results, clusters):
    """
    Save detailed results for a category to files.
    """
    # Create category-specific directory
    import os
    os.makedirs(f"results_{category}", exist_ok=True)
    
    # Save top genes
    top_genes.to_csv(f"results_{category}/top_genes.csv", index=False)
    
    # Save high-fitness pairs
    if fitness_results['high_fitness_pairs']:
        pairs_df = pd.DataFrame(fitness_results['high_fitness_pairs'])
        pairs_df.to_csv(f"results_{category}/high_fitness_pairs.csv", index=False)
    
    # Save cluster assignments
    cluster_df = pd.DataFrame({
        'Gene': data.index,
        'Cluster': clusters,
        'Average_Fitness': data.mean(axis=1)
    })
    cluster_df.to_csv(f"results_{category}/cluster_assignments.csv", index=False)
    
    # Save summary statistics
    summary = {
        'category': category,
        'n_conditions': len(data.columns),
        'n_genes': len(data.index),
        'mean_fitness': data.mean().mean(),
        'std_fitness': data.mean().std(),
        'min_fitness': data.min().min(),
        'max_fitness': data.max().max(),
        'n_high_fitness_pairs': len(fitness_results['high_fitness_pairs'])
    }
    
    summary_df = pd.DataFrame([summary])
    summary_df.to_csv(f"results_{category}/summary_stats.csv", index=False)
    
    print(f"Results saved to results_{category}/ directory")

def create_comparison_plots(results):
    """
    Create comparison plots across categories.
    """
    # Prepare data for comparison
    comparison_data = []
    
    for category, result in results.items():
        data = result['data']
        comparison_data.append({
            'Category': category.replace('_', ' ').title(),
            'Mean_Fitness': data.mean().mean(),
            'Std_Fitness': data.mean().std(),
            'N_Conditions': len(data.columns),
            'N_Genes': len(data.index),
            'N_High_Fitness_Pairs': len(result['fitness_results']['high_fitness_pairs'])
        })
    
    comparison_df = pd.DataFrame(comparison_data)
    
    # Create plots
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # Mean fitness by category
    sns.barplot(data=comparison_df, x='Category', y='Mean_Fitness', ax=axes[0,0])
    axes[0,0].set_title('Mean Fitness by Category')
    axes[0,0].tick_params(axis='x', rotation=45)
    
    # Number of conditions by category
    sns.barplot(data=comparison_df, x='Category', y='N_Conditions', ax=axes[0,1])
    axes[0,1].set_title('Number of Conditions by Category')
    axes[0,1].tick_params(axis='x', rotation=45)
    
    # Fitness variability by category
    sns.barplot(data=comparison_df, x='Category', y='Std_Fitness', ax=axes[1,0])
    axes[1,0].set_title('Fitness Variability (Std) by Category')
    axes[1,0].tick_params(axis='x', rotation=45)
    
    # High-fitness pairs by category
    sns.barplot(data=comparison_df, x='Category', y='N_High_Fitness_Pairs', ax=axes[1,1])
    axes[1,1].set_title('Number of High-Fitness Gene Pairs by Category')
    axes[1,1].tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig('category_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Save comparison data
    comparison_df.to_csv('category_comparison.csv', index=False)
    print("Comparison plots saved as 'category_comparison.png'")
    print("Comparison data saved as 'category_comparison.csv'")

def analyze_gene_networks(results, category='carbon_sources', top_n=50):
    """
    Analyze gene interaction networks within a category.
    """
    if category not in results:
        print(f"Category {category} not found in results")
        return
    
    data = results[category]['data']
    fitness_results = results[category]['fitness_results']
    
    # Get top genes
    top_genes_list = results[category]['top_genes']['Gene'].head(top_n).tolist()
    
    # Filter data to top genes
    top_data = data.loc[top_genes_list]
    
    # Calculate correlation matrix
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(top_data.T).T
    correlation_matrix = np.corrcoef(scaled_data)
    
    # Create network visualization
    fig = go.Figure()
    
    # Add nodes (genes)
    for i, gene in enumerate(top_genes_list):
        fig.add_trace(go.Scatter(
            x=[i], y=[0],
            mode='markers+text',
            text=[gene],
            textposition='top center',
            marker=dict(size=10, color='lightblue'),
            showlegend=False
        ))
    
    # Add edges for high correlations
    threshold = 0.8
    for i in range(len(top_genes_list)):
        for j in range(i+1, len(top_genes_list)):
            if abs(correlation_matrix[i, j]) > threshold:
                fig.add_trace(go.Scatter(
                    x=[i, j], y=[0, 0],
                    mode='lines',
                    line=dict(
                        width=abs(correlation_matrix[i, j]) * 5,
                        color='red' if correlation_matrix[i, j] > 0 else 'blue'
                    ),
                    showlegend=False
                ))
    
    fig.update_layout(
        title=f'Gene Network - {category.replace("_", " ").title()} (Top {top_n} genes)',
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        height=600
    )
    
    fig.show()
    
    # Save correlation matrix
    corr_df = pd.DataFrame(
        correlation_matrix, 
        index=top_genes_list, 
        columns=top_genes_list
    )
    corr_df.to_csv(f'results_{category}/gene_correlations.csv')
    
    print(f"Gene correlation matrix saved to results_{category}/gene_correlations.csv")

if __name__ == "__main__":
    # Run targeted analysis on key categories
    print("Starting targeted functional analysis...")
    
    # Analyze specific categories
    results, analyzer = analyze_specific_categories([
        'carbon_sources', 
        'nitrogen_sources', 
        'antibiotics',
        'metal_stress'
    ])
    
    # Create comparison plots
    print("\nCreating comparison plots...")
    create_comparison_plots(results)
    
    # Analyze gene networks for carbon sources
    print("\nAnalyzing gene networks for carbon sources...")
    analyze_gene_networks(results, 'carbon_sources', top_n=30)
    
    print("\nAnalysis complete! Check the results_* directories for detailed outputs.") 