import pandas as pd
import numpy as np

# Load the fitness values file
df = pd.read_csv("fit_organism_Keio.tsv", sep="\t")

print("Dataset shape:", df.shape)
print("\nColumn names:")
print(df.columns.tolist())

print("\nFirst few rows:")
print(df.head())

print("\nColumn types:")
print(df.dtypes)

# Look for patterns in column names
columns = df.columns.tolist()
print("\nColumns containing 'set':")
set_columns = [col for col in columns if 'set' in col.lower()]
print(set_columns)

print("\nColumns containing 'media' or 'condition':")
media_columns = [col for col in columns if any(word in col.lower() for word in ['media', 'condition', 'medium'])]
print(media_columns)



print("\nColumns containing 'sugars' or 'alcohol':")
sugar_columns = [col for col in columns if any(word in col.lower() for word in ['ol ', 'ose '])]
print(sugar_columns)

print("\nColumns containing 'pH' or 'acid':")
ph_columns = [col for col in columns if any(word in col.lower() for word in ['ph', 'acid', 'alkaline'])]
print(ph_columns)

print("\nColumns containing 'carbon' or 'nitrogen':")
nutrient_columns = [col for col in columns if any(word in col.lower() for word in ['carbon', 'nitrogen', 'nutrient', 'nitr'])]
print(nutrient_columns)

# Save column analysis to file
with open('column_analysis.txt', 'w') as f:
    f.write("Column Analysis for fit_organism_Keio.tsv\n")
    f.write("="*50 + "\n\n")
    f.write(f"Total columns: {len(columns)}\n\n")
    
    f.write("All columns:\n")
    for i, col in enumerate(columns):
        f.write(f"{i+1:3d}. {col}\n")
    
    f.write(f"\nSet columns ({len(set_columns)}):\n")
    for col in set_columns:
        f.write(f"- {col}\n")
    
    f.write(f"\nMedia/Condition columns ({len(media_columns)}):\n")
    for col in media_columns:
        f.write(f"- {col}\n")
    

    f.write(f"\nSugar/Alcohol columns ({len(sugar_columns)}):\n")
    for col in sugar_columns:
        f.write(f"- {col}\n")

print("\nColumn analysis saved to 'column_analysis.txt'") 