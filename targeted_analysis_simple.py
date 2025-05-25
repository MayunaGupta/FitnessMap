import pandas as pd
import numpy as np
from functional_clustering import FunctionalGeneClusterAnalyzer, categorize_columns
import os

def analyze_specific_categories(categories_to_analyze=['carbon_sources', 'nitrogen_sources', 'antibiotics']):
    """
    Perform detailed analysis on specific functional categories.
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
    summary_data = []
    
    for category in categories_to_analyze:
        if category in analyzer.category_data and len(analyzer.category_data[category].columns) > 0:
            print(f"\n{'='*60}")
            print(f"DETAILED ANALYSIS: {category.upper().replace('_', ' ')}")
            print(f"{'='*60}")
            
            data = analyzer.category_data[category]
            print(f"Number of conditions: {len(data.columns)}")
            print(f"Number of genes: {len(data.index)}")
            
            # Basic statistics
            mean_fitness = data.mean().mean()
            std_fitness = data.mean().std()
            min_fitness = data.min().min()
            max_fitness = data.max().max()
            
            print(f"\nFitness Statistics:")
            print(f"  Mean fitness: {mean_fitness:.3f}")
            print(f"  Std fitness: {std_fitness:.3f}")
            print(f"  Min fitness: {min_fitness:.3f}")
            print(f"  Max fitness: {max_fitness:.3f}")
            
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
            
            # Store summary data
            summary_data.append({
                'Category': category.replace('_', ' ').title(),
                'N_Conditions': len(data.columns),
                'N_Genes': len(data.index),
                'Mean_Fitness': mean_fitness,
                'Std_Fitness': std_fitness,
                'Min_Fitness': min_fitness,
                'Max_Fitness': max_fitness,
                'N_High_Fitness_Pairs': len(fitness_results['high_fitness_pairs']),
                'Top_Gene': top_genes.iloc[0]['Gene'],
                'Top_Gene_Fitness': top_genes.iloc[0]['Average_Fitness']
            })
            
            results[category] = {
                'data': data,
                'top_genes': top_genes,
                'fitness_results': fitness_results,
                'clusters': clusters
            }
    
    # Create and save summary comparison
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv('category_summary_comparison.csv', index=False)
    
    print(f"\n{'='*60}")
    print("SUMMARY COMPARISON ACROSS CATEGORIES")
    print(f"{'='*60}")
    print(summary_df.to_string(index=False))
    
    return results, analyzer, summary_df

def save_category_results(category, data, top_genes, fitness_results, clusters):
    """
    Save detailed results for a category to files.
    """
    # Create category-specific directory
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

def analyze_gene_correlations(results, category='carbon_sources', top_n=20):
    """
    Analyze gene correlations within a category.
    """
    if category not in results:
        print(f"Category {category} not found in results")
        return
    
    data = results[category]['data']
    
    # Get top genes
    top_genes_list = results[category]['top_genes']['Gene'].head(top_n).tolist()
    
    # Filter data to top genes
    top_data = data.loc[top_genes_list]
    
    # Calculate correlation matrix
    correlation_matrix = top_data.T.corr()
    
    # Save correlation matrix
    correlation_matrix.to_csv(f'results_{category}/gene_correlations.csv')
    
    # Find highly correlated gene pairs
    high_corr_pairs = []
    threshold = 0.8
    
    for i in range(len(top_genes_list)):
        for j in range(i+1, len(top_genes_list)):
            corr = correlation_matrix.iloc[i, j]
            if abs(corr) > threshold:
                high_corr_pairs.append({
                    'gene1': top_genes_list[i],
                    'gene2': top_genes_list[j],
                    'correlation': corr
                })
    
    # Sort by absolute correlation
    high_corr_pairs.sort(key=lambda x: abs(x['correlation']), reverse=True)
    
    print(f"\nHighly correlated gene pairs in {category} (|correlation| > {threshold}):")
    for i, pair in enumerate(high_corr_pairs[:10]):
        print(f"  {i+1:2d}. {pair['gene1']} - {pair['gene2']}: {pair['correlation']:.3f}")
    
    # Save correlation pairs
    if high_corr_pairs:
        corr_pairs_df = pd.DataFrame(high_corr_pairs)
        corr_pairs_df.to_csv(f'results_{category}/high_correlation_pairs.csv', index=False)
    
    print(f"Gene correlation analysis saved to results_{category}/")
    
    return correlation_matrix, high_corr_pairs

def create_functional_insights(results):
    """
    Generate functional insights from the analysis.
    """
    insights = []
    
    for category, result in results.items():
        top_genes = result['top_genes']['Gene'].head(5).tolist()
        high_fitness_pairs = result['fitness_results']['high_fitness_pairs'][:3]
        
        insight = {
            'category': category,
            'key_findings': {
                'top_genes': top_genes,
                'top_pairs': [(p['gene1'], p['gene2'], p['combined_fitness']) 
                             for p in high_fitness_pairs],
                'n_conditions': len(result['data'].columns),
                'mean_fitness': result['data'].mean().mean()
            }
        }
        insights.append(insight)
    
    # Save insights
    with open('functional_insights.txt', 'w') as f:
        f.write("FUNCTIONAL GENE FITNESS ANALYSIS - KEY INSIGHTS\n")
        f.write("=" * 60 + "\n\n")
        
        for insight in insights:
            category = insight['category'].replace('_', ' ').title()
            f.write(f"{category}:\n")
            f.write(f"  Conditions analyzed: {insight['key_findings']['n_conditions']}\n")
            f.write(f"  Mean fitness: {insight['key_findings']['mean_fitness']:.3f}\n")
            f.write(f"  Top genes: {', '.join(insight['key_findings']['top_genes'])}\n")
            f.write(f"  Top gene pairs:\n")
            for gene1, gene2, fitness in insight['key_findings']['top_pairs']:
                f.write(f"    - {gene1} & {gene2} (fitness: {fitness:.3f})\n")
            f.write("\n")
    
    print("Functional insights saved to 'functional_insights.txt'")
    return insights

if __name__ == "__main__":
    # Run targeted analysis on key categories
    print("Starting targeted functional analysis...")
    
    # Analyze specific categories
    results, analyzer, summary_df = analyze_specific_categories([
        'carbon_sources', 
        'nitrogen_sources', 
        'antibiotics',
        'metal_stress'
    ])
    
    # Analyze gene correlations for each category
    print("\nAnalyzing gene correlations...")
    for category in results.keys():
        analyze_gene_correlations(results, category, top_n=15)
    
    # Generate functional insights
    print("\nGenerating functional insights...")
    insights = create_functional_insights(results)
    
    print("\nAnalysis complete! Check the following outputs:")
    print("- results_*/ directories for detailed category-specific results")
    print("- category_summary_comparison.csv for cross-category comparison")
    print("- functional_insights.txt for key findings summary") 