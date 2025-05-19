# -----------------------------------------------------------------------------
# Script: Phylogenetic Tree Construction from ANIm Data
# Author: [Tom M. Raaymakers] # Please update with the correct author name
# Date: 2025-05-19
# Description: This script takes Average Nucleotide Identity (ANI) percentage
#              data, calculates pairwise distances, performs UPGMA
#              (Unweighted Pair Group Method with Arithmetic Mean)
#              clustering, and outputs the resulting phylogenetic tree
#              in Newick (.nwk) format.
# -----------------------------------------------------------------------------

# --- 1. Package Management ---

# Install 'ape' package if it's not already installed
# 'ape' (Analyses of Phylogenetics and Evolution) is used for handling
# phylogenetic trees.
if (!requireNamespace("ape", quietly = TRUE)) {
  install.packages("ape")
}

# Load the 'ape' library
library(ape)

# --- 2. Set Working Directory (User-specific) ---

setwd("X:/ANIm_output")
# Alternatively, ensure your input file is in the current working directory.

# --- 3. Load and Process Data ---

# Define the input filename
input_file <- "ANIm_percentage_identity.tab"

# Load the ANIm data
ani_data <- read.table(input_file, header = TRUE, row.names = 1)

# Convert ANI percentage identity values to a distance matrix
dist_matrix <- as.dist((1 - ani_data)*100)

# --- 4. Perform Hierarchical Clustering ---

# Perform UPGMA clustering.
upgma_tree <- hclust(dist_matrix, method = "average")

# --- 5. Convert to Phylogenetic Tree Object ---

# Convert the 'hclust' object to a 'phylo' object,
# which is used by 'ape' for phylogenetic trees.
phylo_tree <- as.phylo(upgma_tree)

# --- 6. Output Tree ---

# Define the output filename for the Newick tree
output_newick_file <- "pyani_ANIm_upgma_tree.nwk"

# Write the phylogenetic tree to a Newick format file.
# This format is standard for representing phylogenetic trees.
write.tree(phylo_tree, file = output_newick_file)

# --- 7. Confirmation Message (Optional) ---
cat("UPGMA tree construction complete.\n")
cat("Tree saved to:", output_newick_file, "\n")

# End of script
