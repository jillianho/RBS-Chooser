"""
Example usage of RBSChooser

This notebook demonstrates the basic workflow for selecting optimal
ribosomal binding sites using RBSChooser.
"""

# First, install the package:
# pip install -e .

# Then import the main class
from rbschooser import RBSChooser, validate_cds

# ============================================================================
# Step 1: Initialize RBSChooser with genomic data
# ============================================================================

chooser = RBSChooser()

# Load genomic annotations and proteomics data
# This preprocesses the data once and caches it for reuse
RBSChooser.initiate(
    genbank_file="E_coli_MG1655.gb",           # GenBank genome file
    proteomics_file="proteome_data.txt",       # Expression levels
    top_percent=0.05                           # Use top 5% expressed genes
)

print(f"Initialized with {len(RBSChooser.options)} candidate RBS options")

# ============================================================================
# Step 2: Use the chooser to select RBS for your target genes
# ============================================================================

# Define your target coding sequences
target_genes = {
    "my_gene_1": "ATGGCCATTGTAATGGGCCGCTGCAAGGGT",
    "my_gene_2": "ATGATGATGATGATGATGATGA",
    "my_gene_3": "ATGCTAGCTAGCTAG",
}

# You can also validate a CDS before using it
try:
    validate_cds("ATGATGATG")  # Will raise ValueError (too short)
except ValueError as e:
    print(f"Validation error: {e}")

# For each target gene, find the best RBS
for gene_name, cds in target_genes.items():
    try:
        best_rbs = chooser.run(cds)
        print(f"\n{gene_name}:")
        print(f"  Best RBS from: {best_rbs.gene_name}")
        print(f"  RBS sequence: {best_rbs.utr}")
        print(f"  First 6 AAs: {best_rbs.first_six_aas}")
    except ValueError as e:
        print(f"\n{gene_name}: Invalid CDS - {e}")

# ============================================================================
# Step 3: Exclude specific genes if needed
# ============================================================================

target_cds = "ATGGCCATTGTAATGGGCCGCTGCAAGGGT"

# Get the best option
best = chooser.run(target_cds)
print(f"\nFirst choice: {best.gene_name}")

# Exclude it and find the second-best
second_best = chooser.run(target_cds, ignores={best.gene_name})
print(f"Second choice: {second_best.gene_name}")

# You can exclude multiple genes
third_best = chooser.run(target_cds, ignores={best.gene_name, second_best.gene_name})
print(f"Third choice: {third_best.gene_name}")

# ============================================================================
# Step 4: Verify the selection
# ============================================================================

# The selected RBS option contains useful information
best_rbs = chooser.run("ATGATGATGATGATGATGATGATGATGTAA")

print(f"\nSelected RBS details:")
print(f"  Gene name: {best_rbs.gene_name}")
print(f"  UTR (50bp + RBS): {best_rbs.utr}")
print(f"  CDS start: {best_rbs.cds[:30]}...")
print(f"  First 6 amino acids: {best_rbs.first_six_aas}")

# Construct the optimized construct: RBS + target CDS
target_optimized = best_rbs.utr + target_cds
print(f"\nOptimized construct length: {len(target_optimized)} bp")

# ============================================================================
# Step 5: Work with utility functions directly
# ============================================================================

from rbschooser import (
    reverse_complement,
    first_six_amino_acids,
    occlusion_score,
    peptide_distance,
)

# Reverse complement a sequence
seq = "ATGCGATCG"
rev_comp = reverse_complement(seq)
print(f"\nReverse complement of {seq}: {rev_comp}")

# Get the first 6 amino acids of a CDS
cds = "ATGGCCATTGTAATGGGCCGCTGCAAGGGT"
first_six = first_six_amino_acids(cds)
print(f"First 6 AAs of {cds[:10]}...: {first_six}")

# Check structural interference
utr = "AGGAGGUAAAA"
score = occlusion_score(utr)
print(f"Occlusion score for RBS: {score}")

# Compare peptides
pep_dist = peptide_distance("MATVGR", "MATVGS")
print(f"Peptide distance: {pep_dist}")
