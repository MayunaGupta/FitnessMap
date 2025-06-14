Running Alphafold 2 to predict protein-protein interactions for top gene pairs
The following scripts need to be run in order:
1.	250601_Get_protein_sequences_for_gene_pairs.ipynb
2.	250601_AlphaFold2_Protein_Interactions.ipynb
3.	Get_min_PAE_interaction_Alphafold_prediction.ipynb

1. Takes the following inputs:
- 'ecoli_gls_coefficients.csv’, which is obtained from running GLS.py
-'fit_organism_Keio.tsv', which is downloaded from the E.coli page of the Fitness Browser database. The download link is titled ‘Fitness values (tab-delimited)’ (https://fit.genomics.lbl.gov/cgi-bin/org.cgi?orgId=Keio). 
-‘organism_Keio.faa’, which is also downloaded from the E.coli page of Fitness Browser. The download link is titled ‘Protein sequences (fasta)’
-‘uniprotkb_ecoli_AND_model_organism_8333_2025_06_02.tsv’, which provides a list of all the genes and protein sequences in E.coli from UniProt. This is needed for protein sequences which are missing in ‘organism_Keio.faa’. On UniprotKB ("https://www.uniprot.org/uniprotkb?dir=ascend&facets=model_organism%3A83333&query=ecoli&sort=organism_name" search ecoli and download as a tsv, making sure to include the Gene Names and the Sequence columns.
****All these files are available at the shared google drive folder https://drive.google.com/drive/folders/1-CfFsDSErSE719UJ20IgW3lo0sXQXNh1?usp=drive_link 

The output is a data frame containing the protein sequences for each gene for which a protein sequence is known. Note that some genes do not have protein sequences. These may be RNA encoding, pseudogenes, etc. When I ran the script for the top 50 gene pairs, the genes without protein sequences all encoded RNAs.

2. Takes the following inputs:
-"top_interaction_sequences.csv"
-You need to edit the indexes in the following code block 
for gene_pair in AF_inputs[3:6]:   #predicts protein interaction for the gene pairs of index 3 through 5. Edit this index for the gene pair(s) you want Alphafold to make a prediction for
     runAF(gene_pair[0],gene_pair[1])

IMPORTANT: This is made to run as a google colab. The Alphafold prediction output files for each gene pair will be placed in a zip file and downloaded to your google drive. After obtaining the zip file(s), please upload them to the following shared folder titled 'Alphafold Outputs': https://drive.google.com/drive/folders/1Gz7YlvsNBLRM48aVSyd6X7ODeR3TPklX?usp=drive_link

It takes a very long time to run each prediction, which may be because I only have access to the T4 gpus instead of the A100 gpus. You can run multiple colab pages at once to do more predictions. Let me know if you have any ideas of how to speed things up.

3. Takes the following inputs:
- "top_interaction_sequences.csv"
- "least_interaction_sequences.csv"
- Path to the Alphafold_2_outputs folder that contains all the .zip folders with outputs, each .zip folder has a gene pair’s protein interaction prediction
- Path to an Alphafold 3 outputs folder, if predictions obtained from the AlphaFold server (https://alphafoldserver.com/)
-Min PAE of interaction less than 1.5 is good evidence for a protein-protein interaction. Anything less than 10 might be a possible interaction.

This script will add the Min PAE to the corresponding gene pair row in the data frame. 
