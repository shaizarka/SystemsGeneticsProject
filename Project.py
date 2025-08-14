# --- Summary Section ---
# Data used:
# - QTLs: Data/top_10_snps_QTLs.csv (from HW2)
# - eQTLs: Data/best_eqtls.csv (from HW3)
#
# In HW2, QTL analysis was performed to identify SNPs associated with phenotypic traits (QTLs).
# In HW3, eQTL analysis was performed to identify SNPs associated with gene expression (eQTLs).
# Here, we combine these results to identify SNPs that are both QTLs and eQTLs, and test for causality.

"""
possible relations for causality:

M1: 
L ---> R ---> C

M2:
L ---> C ---> R

M3:
L ---> C
 |
  ---> R

Where L is locus (SNP), R is the gene, and C is the complex trait (phenotype).

"""

import pandas as pd
from utils import get_gene_positions
import math
import scipy.stats
import numpy as np
from itertools import combinations
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings('ignore')
from tqdm import tqdm

SIGNIFICANCE_THRESHOLD = 0.05


def compute_P_Xi_given_Yi(correlation_coef, myu_x, myu_y, sigma_x, sigma_y, Xi, Yi, i):
    """Compute P(Xi|Yi) assuming bivariate normal distribution"""
    if abs(correlation_coef) >= 1.0 or sigma_x == 0 or sigma_y == 0:
        return 1e-10  # Avoid division by zero or invalid correlation
    
    P_Xi_given_Yi = (1/(math.sqrt(2*math.pi*(sigma_x**2)*(1-correlation_coef**2)))) * \
                    math.exp((-(Xi - myu_x - correlation_coef*(sigma_x/sigma_y)*(Yi-myu_y))**2)/ (2*(sigma_x**2)*(1-correlation_coef**2)))
    return max(P_Xi_given_Yi, 1e-10)  # Avoid zero probabilities

def compute_P_X_given_Y(X, Y):
    """Compute P(X|Y) for all elements in X and Y"""
    P_X_given_Y = []
    n = len(X)

    # Flatten values before computing correlation
    correlation_coef = scipy.stats.pearsonr(X.values.flatten(), Y.values.flatten())[0]
    if np.isnan(correlation_coef):
        correlation_coef = 0.0

    # Fix: Use flattened arrays to get scalars
    sigma_x = np.std(X.values.flatten())
    sigma_y = np.std(Y.values.flatten())
    myu_x = np.mean(X.values.flatten())
    myu_y = np.mean(Y.values.flatten())

    # Convert to list of scalars
    X = X.values.flatten().tolist()
    Y = Y.values.flatten().tolist()

    for i in range(n):
        Xi = X[i]
        Yi = Y[i]
        P_Xi_given_Yi = compute_P_Xi_given_Yi(correlation_coef, myu_x, myu_y, sigma_x, sigma_y, Xi, Yi, i)
        P_X_given_Y.append(P_Xi_given_Yi)

    return P_X_given_Y

def compute_P_Li(Li):
    """
    Li is a single SNP genotype value (0, 1, or 2)
    Compute P(Li) assuming mendelian inheritance
    """
    if Li == 0:
        return 0.25  # Homozygous for reference allele
    elif Li == 1:
        return 0.5   # Heterozygous
    elif Li == 2:
        return 0.25  # Homozygous for alternate allele
    else:
        raise ValueError(f"Invalid genotype value: {Li}. Expected 0, 1, or 2.")


def m1_log_likelihood(L, C, R):
    """
    Calculate likelihood for M1: L ---> R ---> C
    P(L, C, R) = P(L) * P(R|L) * P(C|R)
    """
    P_R_given_L = compute_P_X_given_Y(R, L)
    P_C_given_R = compute_P_X_given_Y(C, R)
    
    log_likelihood = 0
    n = len(L)
    L = L.values.flatten()  # Ensure L is a flat array

    for i in range(n):
        # Use log probabilities to avoid underflow
        log_likelihood += np.log(compute_P_Li(L[i])) + np.log(P_R_given_L[i]) + np.log(P_C_given_R[i])

    return log_likelihood

def m2_log_likelihood(L, C, R):
    """
    Calculate likelihood for M2: L ---> C ---> R
    P(L, C, R) = P(L) * P(C|L) * P(R|C)
    """
    P_C_given_L = compute_P_X_given_Y(C, L)
    P_R_given_C = compute_P_X_given_Y(R, C)
    
    log_likelihood = 0
    n = len(L)
    L = L.values.flatten()  # Ensure L is a flat array
    for i in range(n):
        log_likelihood += np.log(compute_P_Li(L[i])) + np.log(P_C_given_L[i]) + np.log(P_R_given_C[i])
    
    return log_likelihood

def m3_log_likelihood(L, C, R):
    """
    Calculate likelihood for M3: L ---> C, L ---> R
    P(L, C, R) = P(L) * P(C|L) * P(R|L)
    """
    P_C_given_L = compute_P_X_given_Y(C, L)
    P_R_given_L = compute_P_X_given_Y(R, L)

    log_likelihood = 0
    n = len(L)
    L = L.values.flatten()  # Ensure L is a flat array
    for i in range(n):
        log_likelihood += np.log(compute_P_Li(L[i])) + np.log(P_C_given_L[i]) + np.log(P_R_given_L[i])

    return log_likelihood

def get_data_for_triplet(snp, gene, phenotype):
    """
    Extract the actual data vectors for a given SNP, gene, and phenotype
    """
    # Get genotype data for the SNP (rows of genotype_df are SNPs)
    if snp in genotype_df["Locus"].values:
        L = genotype_df.loc[genotype_df["Locus"] == snp].iloc[:, 3:]
    else:
        print(f"Warning: QTL SNP {snp} not found in genotype data")
        return None, None, None
    
    # Get expression data for the gene
    if gene in expression_df.columns:
        R = expression_df[[gene]].T
    else:
        print(f"Warning: Gene {gene} not found in expression data")
        return None, None, None
    
    # Get phenotype data
    if phenotype in phenotype_df['ID_FOR_CHECK'].values:
        C = phenotype_df.loc[phenotype_df['ID_FOR_CHECK'] == phenotype].iloc[:, 7:]
    else:
        print(f"Warning: Phenotype {phenotype} not found in phenotype data")
        return None, None, None
    
    # remove any rows with NaN values
    L = L.dropna(axis=1, how='any')
    R = R.dropna(axis=1, how='any')
    C = C.dropna(axis=1, how='any')

    # Ensure all vectors have the same length
    min_len = min(len(L), len(R), len(C))
    L = L[:min_len]
    R = R[:min_len]
    C = C[:min_len]

    l_cols = set(L.columns)
    r_cols = set(R.columns)
    c_cols = set(C.columns)

    # Find common columns
    common_cols = list(l_cols & r_cols & c_cols)

    # Sort them if you want consistent order
    common_cols.sort()

    # Subset each to common columns
    L_common = L[common_cols]
    R_common = R[common_cols]
    C_common = C[common_cols]

    
    return L_common, R_common, C_common

def test_causality(qtl_snp, gene, phenotype):
    """
    Test causality for a given QTL SNP, gene, and phenotype triplet
    """
    L, R, C = get_data_for_triplet(qtl_snp, gene, phenotype)
    
    if L is None or R is None or C is None:
        return None
    
    # Calculate likelihoods for each model
    ll_m1 = m1_log_likelihood(L, C, R)
    ll_m2 = m2_log_likelihood(L, C, R)
    ll_m3 = m3_log_likelihood(L, C, R)
    
    # Find the best model (highest likelihood)
    log_likelihoods = {'M1': ll_m1, 'M2': ll_m2, 'M3': ll_m3}
    best_model = max(log_likelihoods, key=log_likelihoods.get)
    
    return {
        'QTL_SNP': qtl_snp,
        'Gene': gene,
        'Phenotype': phenotype,
        'M1_LogLikelihood': ll_m1,
        'M2_LogLikelihood': ll_m2,
        'M3_LogLikelihood': ll_m3,
        'Best_Model': best_model,
        'Best_LogLikelihood': log_likelihoods[best_model]
    }

def permutation_test(L, R, C, test_model, n_permutations=10_000):
    """
    Perform permutation test to assess statistical significance of causality models

    test_model: str, one of 'M1', 'M2', 'M3'
    
    H0: test_model is not closer to true model than the others
    H1: test_model is closer to true model than the others
    """

    if not test_model in ['M1', 'M2', 'M3']:
        raise ValueError("test_model must be one of 'M1', 'M2', 'M3'")

    # Calculate observed likelihoods
    obs_ll_m1 = m1_log_likelihood(L, C, R)
    obs_ll_m2 = m2_log_likelihood(L, C, R)
    obs_ll_m3 = m3_log_likelihood(L, C, R)
    

    lratios_perm = []
    
    if test_model == 'M1':
        # M1: L -> R -> C
        lratio_obs = obs_ll_m1 - max(obs_ll_m2, obs_ll_m3)
    elif test_model == 'M2':
        # M2: L -> C -> R
        lratio_obs = obs_ll_m2 - max(obs_ll_m1, obs_ll_m3)
    else:
        # M3: L -> C, L -> R
        lratio_obs = obs_ll_m3 - max(obs_ll_m1, obs_ll_m2)
    

    for i in range(n_permutations):
        # Permute the phenotype values (breaking causal relationships)
        if test_model == 'M1':
            C = pd.DataFrame(np.random.permutation(C.values.T).T, index=C.index, columns=C.columns)
        elif test_model == 'M2':
            R = pd.DataFrame(np.random.permutation(R.values.T).T, index=R.index, columns=R.columns)
        else:
            C = pd.DataFrame(np.random.permutation(C.values.T).T, index=C.index, columns=C.columns)

        # Calculate likelihoods with permuted data
        ll_m1_perm = m1_log_likelihood(L, C, R)
        ll_m2_perm = m2_log_likelihood(L, C, R)
        ll_m3_perm = m3_log_likelihood(L, C, R)
        
        if test_model == 'M1':
            lratios_perm.append(ll_m1_perm - max(ll_m2_perm, ll_m3_perm))
        elif test_model == 'M2':
            lratios_perm.append(ll_m2_perm - max(ll_m1_perm, ll_m3_perm))
        else:
            lratios_perm.append(ll_m3_perm - max(ll_m1_perm, ll_m2_perm))

    # Calculate p-values
    p = (np.sum(np.array(lratios_perm) >= lratio_obs)) / (n_permutations)


    return {
        'Original_M1': obs_ll_m1,
        'Original_M2': obs_ll_m2,
        'Original_M3': obs_ll_m3,
        'P_value': p,
        }

# Load QTL and eQTL data
qtl_df = pd.read_csv('Data/top_10_snps_QTLs.csv')
eqtl_df = pd.read_csv('Data/best_eqtls.csv')
phenotype_df = pd.read_csv('Data/phenotypes.csv')
genotype_df = pd.read_csv('Data/genotypes_organized.csv')
expression_df = pd.read_csv('Data/expression_df.csv', index_col=0)  # Assuming expression data exists

phenotype_df = phenotype_df.loc[phenotype_df['ID_FOR_CHECK']=="1446"]
# Drop any genotype columns (strands) that contain 'U' or 'H'
# If metadata columns exist (e.g., first 3), isolate data first
geno_meta = genotype_df.iloc[:, :3]    # adjust if metadata columns differ
geno_data = genotype_df.iloc[:, 3:]

# Find columns that contain any 'U' or 'H'
bad_strands = geno_data.columns[geno_data.isin(['U']).any(axis=0)]

# Drop those columns
geno_data_clean = geno_data.drop(columns=bad_strands)

# Replace B/D with 0/1
geno_data_clean = geno_data_clean.replace({'B': 0, 'H': 1, 'D': 2})

# Reassemble the full genotype dataframe (if you want to preserve metadata)
genotype_df = pd.concat([geno_meta, geno_data_clean], axis=1)

# Find QTL/eQTL pairs that are nearby (within 1Mb)
nearby_pairs = []
window = 1_000_000  # 1Mb
for _, qtl_row in qtl_df.iterrows():
    qtl_snp = qtl_row['SNP']
    qtl_pos = get_gene_positions(qtl_snp, is_snp=True)
    if pd.isnull(qtl_pos):
        continue
    for _, eqtl_row in eqtl_df.iterrows():
        eqtl_snp = eqtl_row['SNP']
        eqtl_pos = get_gene_positions(eqtl_snp, is_snp=True)
        if pd.isnull(eqtl_pos):
            continue
        # Check if on same chromosome and within window
        if abs(qtl_pos - eqtl_pos) <= window:
            nearby_pairs.append({
                'QTL_SNP': qtl_snp,
                'QTL_Fstat': qtl_row['F-stat'],
                'QTL_P': qtl_row['P-value'],
                'eQTL_SNP': eqtl_snp,
                'Gene': eqtl_row['Gene'],
                'eQTL_Fstat': eqtl_row['F-stat'],
                'eQTL_P': eqtl_row['P-value'],
                'Distance': abs(qtl_pos - eqtl_pos),
                'Tissue': eqtl_row['Tissue'],
            })

print(f"Found {len(nearby_pairs)} QTL/eQTL pairs that are nearby (within 1Mb).")
for pair in nearby_pairs:
    print(pair)

# pairs of gene and phenotype (phenotype is always 1446)
gene_phenotype_pairs = []
phenotype = "1446"
for pair in nearby_pairs:
    gene_tissue = f"{pair['Gene']}_{pair['Tissue']}"
    gene_phenotype_pairs.append((gene_tissue, phenotype))
# Remove duplicates
gene_phenotype_pairs = list(set(gene_phenotype_pairs))

print("\nGene/Phenotype pairs (phenotype is always 1446):")
for gene, pheno in gene_phenotype_pairs:
    print(f"Gene: {gene}, Phenotype: {pheno}")


# Apply causality test to all nearby QTL/eQTL pairs
print("\n=== CAUSALITY ANALYSIS RESULTS ===")
causality_results = []

for pair in nearby_pairs:
    qtl_snp = pair['QTL_SNP']
    gene = f"{pair['Gene']}_{pair['Tissue']}"
    phenotype = "1446"
    
    result = test_causality(qtl_snp, gene, phenotype)
    if result:
        causality_results.append(result)
        print(f"\nTriplet: {qtl_snp} - {gene} - {phenotype}")
        print(f"M1 (L->R->C): {result['M1_LogLikelihood']:.2f}")
        print(f"M2 (L->C->R): {result['M2_LogLikelihood']:.2f}")
        print(f"M3 (L->C,R): {result['M3_LogLikelihood']:.2f}")
        print(f"Best Model: {result['Best_Model']}")

# Select 10 triplets for detailed permutation testing
print("\n=== PERMUTATION TESTING ON SELECTED TRIPLETS ===")
print("Selection criteria: Top 10 triplets with highest likelihood differences between models")

# Sort by likelihood difference to select most interesting cases
for result in causality_results:
    likelihoods = [result['M1_LogLikelihood'], result['M2_LogLikelihood'], result['M3_LogLikelihood']]
    result['Likelihood_Range'] = max(likelihoods) - min(likelihoods)

# Sort by likelihood range (most discriminative cases)
causality_results_sorted = sorted(causality_results, key=lambda x: x['Likelihood_Range'], reverse=True)

# make sure triplets in causality_results_sorted are unique
seen = set()
unique_results = []
for res in causality_results_sorted:
    key = f"{res['QTL_SNP']}_{res['Gene']}_{res['Phenotype']}"
    if key not in seen:
        seen.add(key)
        unique_results.append(res)

selected_triplets = unique_results[:min(10, len(unique_results))]

print(f"Running permutation tests on {len(selected_triplets)} triplets...")

permutation_results = []
min_common_stains = float('inf')
for i, result in tqdm(enumerate(selected_triplets), total=len(selected_triplets)):
    qtl_snp = result['QTL_SNP']
    gene = result['Gene']
    phenotype = result['Phenotype']
    
    L, R, C = get_data_for_triplet(qtl_snp, gene, phenotype)
    if all(x is not None for x in (L, R, C)):
        min_common_stains = min(min_common_stains, len(L.columns))

        # Perform permutation test for the best model
        perm_result = permutation_test(L, R, C, result['Best_Model'], n_permutations=10_000)
        perm_result.update({
            'QTL_SNP': qtl_snp,
            'Gene': gene,
            'Phenotype': phenotype,
            'Best_Model': result['Best_Model']
        })
        permutation_results.append(perm_result)

perm_results_df = pd.DataFrame(permutation_results)

print(f"\nMinimum common stains across selected triplets: {min_common_stains}")

if min_common_stains < 7:
    print("Warning: Some triplets have fewer than 7 common stains, which may affect statistical power.")

# FDR correction
fdr_results = multipletests(perm_results_df['P_value'], method='fdr_bh')
perm_results_df['FDR_P_value'] = fdr_results[1]

for _, perm_result in perm_results_df.iterrows():
    print(f"\nPermutation test: {perm_result['QTL_SNP']} - {perm_result['Gene']} - {perm_result['Phenotype']}")
    print(f"Best Model: {perm_result['Best_Model']}")
    print(f"P-value: {perm_result['P_value']:.6f}")
    print(f"FDR corrected P-value: {perm_result['FDR_P_value']:.6f}")
    print(f"Is Significant (p < {SIGNIFICANCE_THRESHOLD}):  {'Yes' if perm_result['FDR_P_value'] < SIGNIFICANCE_THRESHOLD else 'No'}")

# Summary of results
print("\n=== SUMMARY OF PREDICTED RELATIONS ===")
model_counts = {'M1': 0, 'M2': 0, 'M3': 0}
for result in causality_results:
    model_counts[result['Best_Model']] += 1

print(f"Total triplets analyzed: {len(causality_results)}")
print(f"M1 (SNP → Gene → Phenotype): {model_counts['M1']} cases")
print(f"M2 (SNP → Phenotype → Gene): {model_counts['M2']} cases") 
print(f"M3 (SNP → Gene, SNP → Phenotype): {model_counts['M3']} cases")

# Significant results from permutation testing
if permutation_results:
    print("\n=== STATISTICALLY SIGNIFICANT CAUSALITY RESULTS ===")
    
    for perm_result in permutation_results:
        if perm_result['P_value'] < SIGNIFICANCE_THRESHOLD:
            print(f"\n{perm_result['QTL_SNP']} - {perm_result['Gene']} - {perm_result['Phenotype']}")
            print(f"Best model: {perm_result['Best_Model']}")

print("\n=== DESIGN OF PERMUTATION TEST ===")
print("""
Permutation Test Design:
1. For each triplet (L, R, C), we calculate the original likelihoods for all three models
2. We then permute either the phenotype values (C) or the expression values (R) 1000 times, breaking any causal relationships
3. For each permutation, we recalculate the likelihoods for all models
4. The p-value is the fraction of permutations where the permuted likelihood >= original likelihood
5. This tests the null hypothesis that there is no causal relationship

Choice of 10 triplets:
- Selected triplets with the highest likelihood range (max - min across models)
- This ensures we test the most discriminative cases where models differ substantially
- These are the most likely to show meaningful causal relationships
""")