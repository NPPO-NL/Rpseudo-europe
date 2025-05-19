# -----------------------------------------------------------------------------
# Script: Minimum Spanning Tree Analysis and Visualization
# Author: [Tom M. Raaymakers]
# Date: 2025-05-14 
# Description: This script performs cgMLST data processing, calculates
#              pairwise distances, constructs a minimum spanning tree (MST),
#              and visualizes it using ggplot2 and ggnetwork. Specifically, 
#              this code generates the MST for Clade I (Rose clade)
# -----------------------------------------------------------------------------

# --- 1. Setup ---

# Install pacman if not already installed, then load it
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

# Load required packages (pacman will install them if missing)
pacman::p_load(ape, igraph, ggnetwork, ggplot2, ggrepel, tidyverse, readr)

# Define file paths (MODIFY THESE PATHS AS NEEDED)
working_dir         <- "X:/chewie/cladeI/cgMLST_allele_call_100/"
allele_calls_file   <- "results_alleles.tsv"
metadata_file       <- "X:/metadata.txt" #as a tsv see 3.

# Set working directory
if (dir.exists(working_dir)) {
  setwd(working_dir)
  cat("Working directory set to:", getwd(), "\n")
} else {
  stop("Working directory does not exist: ", working_dir)
}

# --- 2. Load and Prepare Allele Call Data ---

# Load allele call data
cat("Loading allele call data...\n")
calls_raw <- readr::read_tsv(file.path(allele_calls_file), col_types = cols(.default = "c")) # Read all as character first

# Set the first column as row names
# Ensure the first column is suitable for row names (unique identifiers)
if (ncol(calls_raw) > 0 && !any(duplicated(calls_raw[[1]]))) {
  calls <- calls_raw %>%
    column_to_rownames(var = colnames(calls_raw)[1])
} else {
  stop("First column is not suitable for row names (missing, empty, or has duplicates).")
}
cat("Allele call data loaded. Dimensions:", dim(calls), "\n")

# --- 3. Load Metadata ---

cat("Loading metadata...\n")
# Assuming metadata is tab-separated with a header.
metadata <- readr::read_tsv(metadata_file, col_types = cols())
cat("Metadata loaded. Dimensions:", dim(metadata), "\n")

# --- 4. Match and Rename Node Names based on Metadata ---

cat("Matching and renaming node names...\n")
calls_rownames <- rownames(calls)
match_indices <- match(calls_rownames, metadata$chewie2)

# Check for NAs in match_indices (i.e., rownames not found in metadata$chewie2)
if (any(is.na(match_indices))) {
  warning("Some rownames in 'calls' data were not found in 'metadata$chewie2'. These will have NA new names.")
  # Optional: print which ones were not found
  # print(calls_rownames[is.na(match_indices)])
}

matched_isolates <- metadata$Isolate2[match_indices]

# Assign new rownames, handling potential NAs if necessary (e.g., by keeping original or a placeholder)
if (length(matched_isolates) == nrow(calls)) {
  # Check for duplicates in new names if they must be unique
  if(any(duplicated(na.omit(matched_isolates)))){
    warning("Duplicate names found after matching with metadata. Rownames might not be unique.")
  }
  rownames(calls) <- matched_isolates
} else {
  stop("Mismatch in length between current rownames and new names after matching.")
}
cat("Node names updated.\n")


# --- 5. Clean Allele Call Data ---

cat("Cleaning allele call data...\n")
# Remove "INF" prefix (but keep the number behind it)
calls_cleaned <- calls %>%
  mutate(across(everything(), ~ gsub("INF", "", .)))

# Remove all other non-numeric values (e.g., LNF, PLOT3, NIPH) and convert to integer
calls_cleaned <- calls_cleaned %>%
  mutate(across(everything(), ~ suppressWarnings(as.integer(.))))
cat("Allele call data cleaned.\n")

# --- 6. Calculate Pairwise Distances ---
# Using the cleaned data 'calls_cleaned'
cat("Calculating pairwise genetic distances...\n")
# pairwise.deletion = TRUE will ignore NAs for a given pair of sites
# pairwise.deletion = FALSE will result in NA distance if any site has NA for a pair of sequences
genetic_dist <- ape::dist.gene(calls_cleaned, method = "pairwise", pairwise.deletion = FALSE, variance = FALSE)
df_genetic_dist <- as.matrix(genetic_dist)
cat("Pairwise distances calculated.\n")

# --- 7. Create Minimum Spanning Tree (MST) ---
cat("Creating Minimum Spanning Tree...\n")
mst_object <- ape::mst(genetic_dist)
mst_graph <- igraph::graph.adjacency(mst_object, mode = "undirected")
cat("MST created.\n")

# --- 8. Prepare Data for ggnetwork Plot ---
cat("Preparing data for ggnetwork plot...\n")
# Extract edge list and corresponding distances for labels
edge_list_matrix <- igraph::get.edgelist(mst_graph)
edge_distances <- round(df_genetic_dist[edge_list_matrix], 2) # Adjust precision if needed

edge_labels_df <- tibble(
  iso1 = edge_list_matrix[, 1],
  iso2 = edge_list_matrix[, 2],
  distance_label = edge_distances
)

# Create base ggnetwork object
# layout_with_fr is a common choice; others include layout_nicely, layout_as_star, etc.
gg_mst <- ggnetwork(mst_graph, arrow.gap = 0, layout = igraph::layout_with_fr(mst_graph))

# Augment gg_mst with node/edge information and metadata
gg_mst_processed <- gg_mst %>%
  mutate(
    node_edge = ifelse(x == xend & y == yend, "node", "edge"),
    isolate_name = ifelse(node_edge == "node", name, NA_character_)
  )

# Create a lookup table for node values
node_lookup <- gg_mst_processed %>%
  filter(node_edge == "node") %>%
  select(x, y, node_name_lookup = name)

# Map end coordinates of edges back to node names
gg_mst_processed <- gg_mst_processed %>%
  left_join(node_lookup, by = c("xend" = "x", "yend" = "y")) %>%
  mutate(
    # For an edge, 'name' is one node, 'node_name_lookup' is the other (from xend, yend)
    # We use these to merge edge distances
    name1 = ifelse(node_edge == "edge", name, NA_character_),
    name2 = ifelse(node_edge == "edge", node_name_lookup, NA_character_)
  ) %>%
  select(-node_name_lookup) # Clean up temporary column

# Merge edge distance labels
# An edge can be (A,B) or (B,A). We need to check both combinations for merging.
gg_mst_processed <- gg_mst_processed %>%
  left_join(edge_labels_df, by = c("name1" = "iso1", "name2" = "iso2")) %>%
  left_join(edge_labels_df, by = c("name1" = "iso2", "name2" = "iso1"), suffix = c("", "_rev")) %>%
  mutate(
    distance = coalesce(distance_label, distance_label_rev)
  ) %>%
  select(-distance_label, -distance_label_rev)

# Merge metadata for NODES
# 'isolate_name' for nodes, or 'name' if you are sure it's the primary ID for nodes
gg_mst_processed <- gg_mst_processed %>%
  left_join(metadata, by = c("isolate_name" = "Isolate2")) # Ensure 'Isolate2' is the correct column in metadata

# Set metadata attributes to NA for edges, as they apply to nodes
metadata_columns <- c("family", "Clade", "Clade2", "Country", "number", "Matrix") # Add all relevant metadata columns here
for (col_name in metadata_columns) {
  if (col_name %in% names(gg_mst_processed)) {
    gg_mst_processed <- gg_mst_processed %>%
      mutate(!!sym(col_name) := ifelse(node_edge == "edge", NA, !!sym(col_name)))
  }
}
cat("ggnetwork data prepared.\n")

# --- 9. Define Colors and Factor Levels ---
cat("Defining plot aesthetics...\n")
family_colors <- c(
  "Araceae" = "pink", "Cucurbitaceae" = "#67AB05", "Moraceae" = "#FF5733",
  "Rosaceae" = "#b21c0e", "Solanaceae" = "#7a5e8d", "Surface water" = "#0f5e9c",
  "Unknown" = "black", "Zingiberaceae" = "#CCB647"
)

clade_colors <- c( # Not used in the current plot aesthetics, but defined
  "XIII" = "brown", "XII" = "magenta", "XI" = "cyan", "X" = "purple",
  "IX" = "orange", "VIII" = "green", "VII" = "turquoise", "VI" = "pink",
  "V" = "yellow2", "IV" = "red", "III" = "grey", "II" = "olivedrab", "I" = "blue"
)

country_colors <- c( # Not used in the current plot aesthetics, but defined
  "Venezuela" = "brown", "Netherlands (import)" = "orange", "Thailand" = "green",
  "Peru" = "turquoise", "China" = "red", "Bangladesh" = "olivedrab",
  "Israel" = "blue", "South Africa" = "brown"
)

# Ensure Clade2 is a factor with the correct order for legend (if used)
clade_levels <- c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII")
if ("Clade2" %in% names(gg_mst_processed)) {
  gg_mst_processed$Clade2 <- factor(gg_mst_processed$Clade2, levels = clade_levels)
}

# --- 10. Plot the MST ---
cat("Generating MST plot...\n")
mst_plot <- ggplot(gg_mst_processed, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = "black", alpha = 0.5, curvature = 0.0) + 
  geom_edgetext(
    aes(label = distance),
    color = "black",
    size = 5, 
    # vjust = 0, # Adjust vertical position of edge labels
    fill = "white" # can be set to NA
  ) +
  geom_nodes(aes(color = family), size = 13) + 
  geom_nodetext(aes(label = number), color = "black", size = 5) + 
  scale_color_manual(values = family_colors, name = "Family", na.translate = FALSE) + 
  guides(color = guide_legend(override.aes = list(size = 13))) + # Adjusted legend symbol size
  theme_blank() +
  theme(
    legend.position = c(0.85, 0.85), # Adjust as needed
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.key.size = unit(1, "lines")
    )
cat("MST plot generated.\n")
# print(mst_plot) # Uncomment to display plot in R session

# --- 11. Save the Plot ---
current_date <- format(Sys.Date(), "%Y%m%d")
output_filename <- paste0(current_date, "_minimum_spanning_tree_clade_I.pdf")

cat("Saving plot to:", output_filename, "\n")
ggsave(
  filename = output_filename,
  plot = mst_plot,
  device = "pdf",
  width = 10, 
  height = 8, 
  units = "in"
  # path = getwd() # 'path' is optional if already in the desired working directory
)
cat("Script finished successfully.\n")
