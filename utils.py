
import pandas as pd
snp_positions = pd.read_csv("Data/genotypes_organized.csv")[["Locus", "Build37_position"]]
gene_positions = pd.read_csv("Data/MGI_Coordinates.Build37.rpt.txt", delimiter = "\t")[["marker symbol", "representative genome start"]]

def get_gene_positions(gene_name, is_snp=False):
    """
    Returns a genome position (int) for a SNP or gene.
    Returns None if not found.
    """
    if is_snp:
        gene_mask = snp_positions["Locus"] == gene_name
        if gene_mask.any():
            gene = snp_positions[gene_mask].iloc[0]
            return gene["Build37_position"]
        else:
            return None
    else:
        gene_mask = gene_positions["marker symbol"] == gene_name
        if gene_mask.any():
            gene = gene_positions[gene_mask].iloc[0]
            return gene["representative genome start"]
        else:
            return None


def classify_cis_trans(row, threshold=20_000_000):
    try:
        gene_start_pos = get_gene_positions(row["Gene"], is_snp=False)
        snp_start_pos = get_gene_positions(row["SNP"], is_snp=True)

        if abs(gene_start_pos - snp_start_pos) <= threshold:
            return "cis"
        else:
            return "trans"
    except:
        return "unknown"
