# Plot Histogram to show distribution of Number of associations per SNP
import matplotlib.pyplot as plt
from utils import get_gene_positions
import pandas as pd
import numpy as np
import os

def plot_association_distribution(df):
    """
    Plots a histogram showing the distribution of the number of associations per SNP.
    Meaning - How many times a specific SNP appears in the dataframe as associated with a trait or condition.

    Parameters:
    df (pd.DataFrame): DataFrame containing SNP data with a 'SNP' column.
    """

    # Plotting the histogram
    plt.figure(figsize=(10, 6))
    # Dataframe contains record per eQTL with the associated SNPS

    if 'SNP' not in df.columns:
        raise ValueError("DataFrame must contain a 'SNP' column.")
    
    # Count gene name per SNP
    gene_counts = df['SNP'].value_counts()
    plt.hist(gene_counts, bins=30, color='blue', alpha=0.5, edgecolor='black')
    plt.title('Distribution of Number of genes per eQTL')
    plt.xlabel('Number of genes per eQTL')
    plt.ylabel('Frequency')
    plt.grid(axis='y', alpha=0.75)
    plt.tight_layout()
    plt.show()


def plot_association_count_over_genome(df, fig_ax=None, is_scatter=False, should_show=False, color='blue', label=None):
    """
    Plots the number of associations per genomic region.
    """
    if fig_ax is None:
        fig, ax = plt.subplots(figsize=(12, 6))
    else:
        fig, ax = fig_ax
    ax.set_title('Number of Genes per eQTL across Genomic Region')
    ax.set_xlabel('Genomic Position')
    ax.set_ylabel('Number of Genes per eQTL')

    # Count associations per genomic region
    gene_counts = df['SNP'].value_counts()
    # Extract position values directly from get_gene_positions
    # Note: get_gene_positions returns the position directly, not a dictionary
    positions = []
    for snp in gene_counts.index:
        try:
            pos = get_gene_positions(snp, is_snp=True)
            positions.append(pos)
        except Exception as e:
            print(f"Could not get position for {snp}: {e}")
            positions.append(0)
    
    # Sort by position for better visualization
    sorted_data = sorted(zip(positions, gene_counts.values))
    positions_sorted, counts_sorted = zip(*sorted_data) if sorted_data else ([], [])

    if is_scatter:
        ax.vlines(positions_sorted, ymin=0, ymax=counts_sorted, color=color, linestyle='--', alpha=0.5, label=label)
    else:
        ax.plot(positions_sorted, counts_sorted,
                 'o-', color=color,
                 markersize=4, linewidth=1,
                 label=label)
    # ax.set_xticks(rotation=90)
    ax.grid(axis='y', alpha=0.75)

    ax.legend()
    if should_show:
        # Add legend 

        plt.show()

    return fig, ax
    
# def plot_enhanced_association_count_over_genome(df):
#     """
#     Plots the number of associations per genomic region with chromosome markers.
#     This is an enhanced version that divides the genome into chromosomes and 
#     provides a more informative visualization.
#     """

    
#     # Get all chromosomes from genotypes file
#     genotypes_file = "Data/genotypes_organized.csv"
    
#     if not os.path.exists(genotypes_file):
#         print(f"Error: File {genotypes_file} not found!")
#         print("Please check that the genotypes file exists in the Data directory.")
#         return
        
#     try:
#         genotypes = pd.read_csv(genotypes_file)
        
#         # Verify required columns exist
#         required_columns = ['Locus', 'Chr_Build37', 'Build37_position']
#         for col in required_columns:
#             if col not in genotypes.columns:
#                 print(f"Error: Required column '{col}' not found in genotypes file!")
#                 print(f"Available columns are: {', '.join(genotypes.columns)}")
#                 return
        
#         # Create a mapping of SNPs to chromosomes and positions
#         snp_info = {}
#         for _, row in genotypes.iterrows():
#             snp_info[row['Locus']] = {
#                 'chr': row['Chr_Build37'],
#                 'position': row['Build37_position']
#             }
#     except Exception as e:
#         print(f"Error loading or processing genotypes file: {e}")
#         return
    
#     # Count associations per SNP
#     gene_counts = df['SNP'].value_counts()
    
#     # Prepare data for plotting by chromosome
#     chr_data = {}
#     unknown_snps = []
    
#     for snp, count in gene_counts.items():
#         if snp in snp_info:
#             chromosome = snp_info[snp]['chr']
#             position = snp_info[snp]['position']
#             if chromosome not in chr_data:
#                 chr_data[chromosome] = {'positions': [], 'counts': []}
#             chr_data[chromosome]['positions'].append(position)
#             chr_data[chromosome]['counts'].append(count)
#         else:
#             unknown_snps.append(snp)
    
#     if unknown_snps:
#         print(f"Warning: {len(unknown_snps)} SNPs not found in genotypes file")
    
#     # Plot by chromosome
#     chromosomes = sorted(chr_data.keys())
#     n_chromosomes = len(chromosomes)
    
#     # Set up the figure
#     fig, axes = plt.subplots(n_chromosomes, 1, figsize=(12, 2 * n_chromosomes))
#     fig.suptitle('Number of Associations per Genomic Region by Chromosome', fontsize=16)
    
#     # If there's only one chromosome, wrap axes in a list
#     if n_chromosomes == 1:
#         axes = [axes]
    
#     # Plot each chromosome
#     for i, chromosome in enumerate(chromosomes):
#         ax = axes[i]
#         positions = chr_data[chromosome]['positions']
#         counts = chr_data[chromosome]['counts']
        
#         # Sort by position
#         sorted_data = sorted(zip(positions, counts))
#         if sorted_data:
#             positions_sorted, counts_sorted = zip(*sorted_data)
            
#             # Plot as a line with markers
#             ax.plot(positions_sorted, counts_sorted, 'o-', markersize=4)
#             ax.set_title(f'Chromosome {chromosome}')
#             ax.set_xlabel('Position (bp)')
#             ax.set_ylabel('Associations')
#             ax.grid(True, alpha=0.3)
            
#             # Highlight hotspots (regions with high association counts)
#             threshold = np.percentile(counts_sorted, 90)  # Top 10% are hotspots
#             for pos, count in zip(positions_sorted, counts_sorted):
#                 if count >= threshold:
#                     ax.axvline(x=pos, color='r', linestyle='--', alpha=0.3)
    
#     plt.tight_layout(rect=[0, 0, 1, 0.96])
#     plt.show()
    
#     # Also create a summary plot showing distribution across chromosomes
#     plt.figure(figsize=(10, 6))
    
#     # Count total associations per chromosome
#     chr_counts = {chr_name: sum(chr_data[chr_name]['counts']) for chr_name in chr_data}
#     chr_names = list(chr_counts.keys())
#     counts = list(chr_counts.values())
    
#     # Sort by chromosome number
#     chr_pairs = [(int(chr_name) if chr_name.isdigit() else 999, chr_name, count) 
#                  for chr_name, count in zip(chr_names, counts)]
#     chr_pairs.sort()
    
#     sorted_chr_names = [pair[1] for pair in chr_pairs]
#     sorted_counts = [pair[2] for pair in chr_pairs]
    
#     plt.bar(sorted_chr_names, sorted_counts, color='skyblue')
#     plt.xlabel('Chromosome')
#     plt.ylabel('Number of Associations')
#     plt.title('Total eQTL Associations by Chromosome')
#     plt.xticks(rotation=45)
#     plt.grid(axis='y', alpha=0.3)
#     plt.tight_layout()
#     plt.show()