#### Testing GENIE3 (https://github.com/vahuynh/GENIE3/tree/master) to reconstruct tissue-specific GRNs using an unsupervised ML method
### This will regress each gene’s expression on TFs (restricted to footprint‑supported TF–target pairs) and rank TF importance.
### Keep only edges supported by both (a) co‑expression signal and (b) footprint evidence.
### Use Jaccard‑style module overlap or edge‑overlap statistics to quantify how much GRN structure is shared vs. rewired across tissues/species.
### How does this translate in your HA-HE genes and examples discussed in the paper e.g., MAF-actr1?

# STEPS  
# 1. reads *_{B,E,L,T}.expr_mat_st.tsv files,
# 2. applies log₂TPM+1 thresholding (TPM ≥ 1 in at least 2 samples),
# 3. runs GENIE3 per species–tissue, and
# 4. saves footprint‑filtered TF→target edges.




# ================================================
# GENIE3 Gene Regulatory Network Inference
# ngs101 tutorial style + TF filtering (https://ngs101.com/how-to-build-gene-regulatory-networks-from-rna-seq-data-using-genie3-complete-step-by-step-guide-for-absolute-beginners/)
# ================================================

# Install/load packages
required_packages <- c("GENIE3", "data.table", "dplyr", "doParallel","pheatmap","ggdist","readr")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    if (pkg == "GENIE3") {
      if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install("GENIE3")
    } else {
      install.packages(pkg)
    }
  }
  library(pkg, character.only = TRUE)
}

# Paths
DATA_DIR <- "/Users/tarangmehta/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/TGAC/Technology_Development/ATAC-Seq/ATAC_manuscript/Files/genie3"
OUT_DIR  <- DATA_DIR
TF_FILE  <- "TF_Orthogroup_UnifiedGeneSymbol.tsv"

# Parameters (ngs101 style)
MIN_CPM      <- 1        # CPM equivalent threshold for log2(TPM+1) - Removes genes that are silent or barely detectable across most replicates
MIN_SAMPLES  <- 2        # minimum samples with expression above threshold
VAR_PERCENT  <- 0.20     # keep top 80% variable genes - remove the bottom 20% of genes by variance—those that barely change across samples. These genes are unlikely to reveal meaningful regulatory relationships since effective regulation requires expression changes.
MIN_REPS_VAR <- 3        # minimum replicates for variance filtering - ONLY perform variance filtering if more than 2 samples
N_TREES      <- 1000
N_CORES      <- 4

# Tissue codes
tissue_codes <- c("B", "E", "L", "T")

# -----------------------------------
# Load master TF list (column 1 = Gene IDs)
# -----------------------------------
tf_master <- read.table(file.path(DATA_DIR, TF_FILE), sep = "\t", header = FALSE)
tf_master_genes <- unname(tf_master[, 1])  # Fixed: column 1 for your file
cat("Master TF list:", length(tf_master_genes), "TFs\n") # Master TF list: 877 TFs


# -----------------------------------
# AUTO-DETECT all species-tissue files
# -----------------------------------
all_files <- list.files(DATA_DIR, pattern = ".*_[B|E|L|T]\\.expr_mat_st\\.tsv$", full.names = TRUE)
cat("Found", length(all_files), "files to process:\n")
print(basename(all_files))

# Extract species codes from filenames
species_list <- unique(gsub("_(B|E|L|T)\\.expr_mat_st\\.tsv$", "", basename(all_files)))
cat("\nSpecies detected:", paste(species_list, collapse = ", "), "\n")

links_all <- list()

# -----------------------------------
# Process ALL species-tissues
# -----------------------------------
for (species_id in species_list) {
  cat("\n", paste0(rep("=", 50), collapse=""), "\n")
  cat("PROCESSING SPECIES:", species_id, "\n")
  cat(paste0(rep("=", 50), collapse=""), "\n")
  
  species_links <- list()
  
  for (tc in tissue_codes) {
    fn <- file.path(DATA_DIR, sprintf("%s_%s.expr_mat_st.tsv", species_id, tc))
    
    if (!file.exists(fn)) {
      cat("  Skipping", species_id, tc, "(file not found)\n")
      next
    }
    
    cat("\n  ---", species_id, tc, "---\n")
    
    # Read expression matrix
    expr_mat <- fread(fn, header = TRUE, data.table = FALSE)
    gene_ids <- expr_mat[, 1]
    rownames(expr_mat) <- gene_ids
    expr_mat <- as.matrix(expr_mat[, -1])
    
    cat("    Matrix:", nrow(expr_mat), "genes x", ncol(expr_mat), "samples\n")
    
    # Step 1: CPM filter
    keep_genes <- rowSums(expr_mat >= MIN_CPM) >= MIN_SAMPLES
    expr_filtered <- expr_mat[keep_genes, ]
    cat("    After CPM filter:", nrow(expr_filtered), "genes\n")
    
    n_replicates <- ncol(expr_filtered)
    
    # Step 2: Conditional variance filtering
    if (n_replicates >= MIN_REPS_VAR) {
      gene_vars <- apply(expr_filtered, 1, var, na.rm = TRUE)
      var_threshold <- quantile(gene_vars, VAR_PERCENT, na.rm = TRUE)
      keep_var <- gene_vars >= var_threshold
      expr_final <- expr_filtered[keep_var, ]
      cat("    Variance filter (top", (1-VAR_PERCENT)*100, "%):", nrow(expr_final), "genes\n")
    } else {
      expr_final <- expr_filtered
      cat("    Skipped variance filter (n =", n_replicates, "<", MIN_REPS_VAR, ")\n")
    }
    
    # Step 3: Filter TFs to expressed ones only
    available_tfs <- intersect(rownames(expr_final), tf_master_genes)
    cat("    TFs available:", length(available_tfs), "/", length(tf_master_genes), "\n")
    
    if (length(available_tfs) == 0) {
      cat("    WARNING: No TFs expressed, skipping\n")
      next
    }
    
    tf_genes_tissue <- available_tfs
    
    # Step 4: GENIE3
    set.seed(123 + match(species_id, species_list) * 10 + match(tc, tissue_codes))
    
    weight_mat <- GENIE3(
      exprMatrix  = expr_final,
      regulators  = tf_genes_tissue,
      nTrees      = N_TREES,
      nCores      = N_CORES
    )
    
    # Step 5: Extract links
    links <- getLinkList(weight_mat)
    colnames(links) <- c("Regulator", "Target", "Weight")
    links <- links[!is.na(links$Regulator) & !is.na(links$Target), ]
    links <- links[order(links$Weight, decreasing = TRUE), ]
    
    cat("    Generated", nrow(links), "TF→target links\n")
    
    # Step 6: Save
    out_fn <- file.path(OUT_DIR, sprintf("%s_%s_GENIE3_links.tsv", species_id, tc))
    write.table(links, out_fn, sep = "\t", quote = FALSE, row.names = FALSE)
    
    species_links[[tc]] <- links
    cat("    Saved:", basename(out_fn), "\n")
  }
  
  links_all[[species_id]] <- species_links
}


## SUMMARY OF GENIE3 GRNS OUTPUTS:
# ================================================== 
#   PROCESSING SPECIES: Ab 
# ================================================== 
#   
#   --- Ab B ---
#   Matrix: 14094 genes x 2 samples
# After CPM filter: 7212 genes
# Skipped variance filter (n = 2 < 3 )
# TFs available: 221 / 877 
# Generated 1593631 TF→target links
# Saved: Ab_B_GENIE3_links.tsv 
# 
# --- Ab E ---
#   Matrix: 10254 genes x 2 samples
# After CPM filter: 6912 genes
# Skipped variance filter (n = 2 < 3 )
# TFs available: 211 / 877 
# Generated 1458221 TF→target links
# Saved: Ab_E_GENIE3_links.tsv 
# 
# --- Ab L ---
#   Matrix: 9838 genes x 2 samples
# After CPM filter: 2239 genes
# Skipped variance filter (n = 2 < 3 )
# TFs available: 50 / 877 
# Generated 111900 TF→target links
# Saved: Ab_L_GENIE3_links.tsv 
# 
# --- Ab T ---
#   Matrix: 18420 genes x 2 samples
# After CPM filter: 5516 genes
# Skipped variance filter (n = 2 < 3 )
# TFs available: 150 / 877 
# Generated 827250 TF→target links
# Saved: Ab_T_GENIE3_links.tsv 
# 
# ================================================== 
#   PROCESSING SPECIES: Mz 
# ================================================== 
#   
#   --- Mz B ---
#   Matrix: 16970 genes x 2 samples
# After CPM filter: 4462 genes
# Skipped variance filter (n = 2 < 3 )
# TFs available: 131 / 877 
# Generated 584391 TF→target links
# Saved: Mz_B_GENIE3_links.tsv 
# 
# --- Mz E ---
#   Matrix: 16319 genes x 2 samples
# After CPM filter: 2807 genes
# Skipped variance filter (n = 2 < 3 )
# TFs available: 66 / 877 
# Generated 185196 TF→target links
# Saved: Mz_E_GENIE3_links.tsv 
# 
# --- Mz L ---
#   Matrix: 11038 genes x 2 samples
# After CPM filter: 2701 genes
# Skipped variance filter (n = 2 < 3 )
# TFs available: 61 / 877 
# Generated 164700 TF→target links
# Saved: Mz_L_GENIE3_links.tsv 
# 
# --- Mz T ---
#   Matrix: 17878 genes x 2 samples
# After CPM filter: 5789 genes
# Skipped variance filter (n = 2 < 3 )
# TFs available: 150 / 877 
# Generated 868200 TF→target links
# Saved: Mz_T_GENIE3_links.tsv 
# 
# ================================================== 
#   PROCESSING SPECIES: Nb 
# ================================================== 
#   
#   --- Nb B ---
#   Matrix: 16016 genes x 2 samples
# After CPM filter: 6553 genes
# Skipped variance filter (n = 2 < 3 )
# TFs available: 213 / 877 
# Generated 1395576 TF→target links
# Saved: Nb_B_GENIE3_links.tsv 
# 
# --- Nb E ---
#   Matrix: 12149 genes x 2 samples
# After CPM filter: 2334 genes
# Skipped variance filter (n = 2 < 3 )
# TFs available: 90 / 877 
# Generated 209970 TF→target links
# Saved: Nb_E_GENIE3_links.tsv 
# 
# --- Nb L ---
#   Matrix: 12417 genes x 2 samples
# After CPM filter: 5133 genes
# Skipped variance filter (n = 2 < 3 )
# TFs available: 141 / 877 
# Generated 723612 TF→target links
# Saved: Nb_L_GENIE3_links.tsv 
# Skipping Nb T (file not found)
# 
# ================================================== 
#   PROCESSING SPECIES: On 
# ================================================== 
#   
#   --- On B ---
#   Matrix: 18716 genes x 3 samples
# After CPM filter: 7145 genes
# Variance filter (top 80 %): 5716 genes
# TFs available: 180 / 877 
# Generated 1028661 TF→target links
# Saved: On_B_GENIE3_links.tsv 
# 
# --- On E ---
#   Matrix: 19307 genes x 3 samples
# After CPM filter: 4766 genes
# Variance filter (top 80 %): 3813 genes
# TFs available: 92 / 877 
# Generated 350699 TF→target links
# Saved: On_E_GENIE3_links.tsv 
# 
# --- On L ---
#   Matrix: 15675 genes x 3 samples
# After CPM filter: 6299 genes
# Variance filter (top 80 %): 5039 genes
# TFs available: 137 / 877 
# Generated 690173 TF→target links
# Saved: On_L_GENIE3_links.tsv 
# 
# --- On T ---
#   Matrix: 22149 genes x 2 samples
# After CPM filter: 17586 genes
# Skipped variance filter (n = 2 < 3 )
# TFs available: 661 / 877 
# Generated 11623685 TF→target links
# Saved: On_T_GENIE3_links.tsv 
# 
# ================================================== 
#   PROCESSING SPECIES: Pn 
# ================================================== 
#   
#   --- Pn B ---
#   Matrix: 17964 genes x 4 samples
# After CPM filter: 9530 genes
# Variance filter (top 80 %): 7624 genes
# TFs available: 269 / 877 
# Generated 2050572 TF→target links
# Saved: Pn_B_GENIE3_links.tsv 
# 
# --- Pn E ---
#   Matrix: 15807 genes x 3 samples
# After CPM filter: 2876 genes
# Variance filter (top 80 %): 2301 genes
# TFs available: 51 / 877 
# Generated 117294 TF→target links
# Saved: Pn_E_GENIE3_links.tsv 
# 
# --- Pn L ---
#   Matrix: 12268 genes x 4 samples
# After CPM filter: 3763 genes
# Variance filter (top 80 %): 3010 genes
# TFs available: 69 / 877 
# Generated 207621 TF→target links
# Saved: Pn_L_GENIE3_links.tsv 
# 
# --- Pn T ---
#   Matrix: 19166 genes x 4 samples
# After CPM filter: 9165 genes
# Variance filter (top 80 %): 7332 genes
# TFs available: 220 / 877 
# Generated 1612785 TF→target links
# Saved: Pn_T_GENIE3_links.tsv


# ================================================
# Complete GENIE3 Comparative GRN Analysis - GENE EXPRESSION ONLY (no footprint info)
# - Reads *_GENIE3_links.tsv
# - Keeps the top 100,000 edges per network
# - Computes tissue-specific Jaccard matrices
# - Plots heatmaps
# - Extracts top rewired edges
# - Maps Orthogroup IDs to Unified_GeneSymbol
# - Builds publication-style multiplots
# ================================================

# -----------------------------
# Packages
# -----------------------------
required_packages <- c(
  "data.table", "dplyr", "tidyr", "ggplot2", "patchwork",
  "pheatmap", "stringr", "purrr", "readr", "png", "tibble"
)

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

# -----------------------------
# Paths and parameters
# -----------------------------
DATA_DIR <- "/Users/tarangmehta/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/TGAC/Technology_Development/ATAC-Seq/ATAC_manuscript/Files/genie3"
OUT_DIR  <- file.path(DATA_DIR, "comparative_GRN")
ORTHOLOGS_FILE <- "/Users/tarangmehta/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/TGAC/Technology_Development/ATAC-Seq/ATAC_manuscript/Files/OrthologousGroups_ENS_GeneSymbols_6cich.txt"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Compare only same tissues across species
tissue_codes <- c("B", "E", "L", "T")

# Tissue names
tissue_names <- c(
  "B" = "Brain",
  "E" = "Retina",
  "L" = "Liver",
  "T" = "Testis"
)

# Keep top N edges from each network
TOP_N_EDGES <- 12000

# Number of top rewired edges to report
TOP_REWIRED <- 100

# Optional minimum number of species an edge must appear in to be called shared
MIN_SHARED_SPECIES <- 2

# -----------------------------
# Helper functions
# -----------------------------
get_edge_id <- function(df) {
  paste(df$Regulator, df$Target, sep = "->")
}

jaccard_index <- function(edges1, edges2) {
  inter <- length(intersect(edges1, edges2))
  union <- length(union(edges1, edges2))
  if (union == 0) return(NA_real_)
  inter / union
}

safe_read_links <- function(file) {
  df <- fread(file, header = TRUE, data.table = FALSE)
  colnames(df)[1:3] <- c("Regulator", "Target", "Weight")
  
  df <- df %>%
    filter(!is.na(Regulator), !is.na(Target), !is.na(Weight)) %>%
    mutate(
      Regulator = as.character(Regulator),
      Target    = as.character(Target),
      Weight    = as.numeric(Weight)
    ) %>%
    arrange(desc(Weight)) %>%
    mutate(
      Rank = row_number(),
      Edge = paste(Regulator, Target, sep = "->")
    )
  
  df
}

map_oma_to_symbol <- function(gene_col, orthologs_dt) {
  out <- data.table(
    Original_Gene = as.character(gene_col),
    Unified_GeneSymbol = NA_character_,
    Mapping_Success = FALSE
  )
  
  idx <- match(out$Original_Gene, orthologs_dt$Orthogroup)
  hit <- !is.na(idx)
  
  out$Unified_GeneSymbol[hit] <- orthologs_dt$Unified_GeneSymbol[idx[hit]]
  out$Mapping_Success[hit] <- !is.na(out$Unified_GeneSymbol[hit]) & out$Unified_GeneSymbol[hit] != ""
  
  out
}

safe_scale <- function(x) {
  if (all(is.na(x))) return(rep(0, length(x)))
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  as.numeric(scale(x))
}

# -----------------------------
# Find all GENIE3 link files
# -----------------------------
all_files <- list.files(
  DATA_DIR,
  pattern = "^[A-Za-z]+_[BELT]_GENIE3_links\\.tsv$",
  full.names = TRUE
)

if (length(all_files) == 0) {
  stop("No *_GENIE3_links.tsv files found in DATA_DIR")
}

file_table <- tibble(
  file = all_files,
  base = basename(all_files)
) %>%
  mutate(
    Species = str_extract(base, "^[A-Za-z]+"),
    Tissue  = str_extract(base, "(?<=_)[BELT](?=_GENIE3_links\\.tsv$)")
  ) %>%
  arrange(Tissue, Species)

write.table(
  file_table,
  file.path(OUT_DIR, "detected_GENIE3_link_files.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# -----------------------------
# Read all networks
# -----------------------------
network_list <- list()

for (i in seq_len(nrow(file_table))) {
  f <- file_table$file[i]
  sp <- file_table$Species[i]
  ts <- file_table$Tissue[i]
  
  cat("Reading:", basename(f), "\n")
  df <- safe_read_links(f)
  df_top <- head(df, min(TOP_N_EDGES, nrow(df)))
  
  network_list[[paste(sp, ts, sep = "_")]] <- list(
    species = sp,
    tissue = ts,
    full = df,
    top = df_top
  )
}

# -----------------------------
# Tissue-specific Jaccard matrices
# -----------------------------
jaccard_summary <- list()

for (tissue in tissue_codes) {
  tissue_nets <- network_list[purrr::map_chr(network_list, "tissue") == tissue]
  
  if (length(tissue_nets) < 2) {
    cat("Skipping tissue", tissue, "- fewer than 2 species available\n")
    next
  }
  
  species_names <- purrr::map_chr(tissue_nets, "species")
  edge_sets <- purrr::map(tissue_nets, ~ .x$top$Edge)
  
  mat <- matrix(
    NA_real_,
    nrow = length(species_names),
    ncol = length(species_names),
    dimnames = list(species_names, species_names)
  )
  
  for (i in seq_along(species_names)) {
    for (j in seq_along(species_names)) {
      mat[i, j] <- jaccard_index(edge_sets[[i]], edge_sets[[j]])
    }
  }
  
  jaccard_summary[[tissue]] <- as.data.frame(mat)
  
  write.table(
    mat,
    file.path(OUT_DIR, paste0("Jaccard_", tissue, "_species_matrix.tsv")),
    sep = "\t", quote = FALSE, col.names = NA
  )
  
  pheatmap::pheatmap(
    mat,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    display_numbers = TRUE,
    number_format = "%.2f",
    color = colorRampPalette(c("white", "#fee08b", "#d73027"))(100),
    border_color = NA,
    main = paste("Jaccard similarity of top", TOP_N_EDGES, "GENIE3 edges -", tissue_names[[tissue]]),
    filename = file.path(OUT_DIR, paste0("Jaccard_", tissue, "_heatmap.png")),
    width = 7,
    height = 6
  )
  
  long_df <- as.data.frame(as.table(mat)) %>%
    rename(Species1 = Var1, Species2 = Var2, Jaccard = Freq) %>%
    arrange(desc(Jaccard))
  
  write.table(
    long_df,
    file.path(OUT_DIR, paste0("Jaccard_", tissue, "_long.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

# -----------------------------
# Distance / MDS plots per tissue
# -----------------------------
for (tissue in names(jaccard_summary)) {
  mat <- as.matrix(jaccard_summary[[tissue]])
  dist_mat <- 1 - mat
  diag(dist_mat) <- 0
  
  mds <- cmdscale(as.dist(dist_mat), k = 2, eig = TRUE)
  mds_df <- data.frame(
    Species = rownames(mat),
    Dim1 = mds$points[, 1],
    Dim2 = mds$points[, 2]
  )
  
  write.table(
    mds_df,
    file.path(OUT_DIR, paste0("MDS_", tissue, "_species_positions.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  p <- ggplot(mds_df, aes(Dim1, Dim2, label = Species)) +
    geom_point(size = 4, color = "#2c7fb8") +
    geom_text(vjust = -0.8, size = 4) +
    theme_bw(base_size = 12) +
    labs(
      title = paste("Species distance in GRN space -", tissue_names[[tissue]]),
      subtitle = "Distance = 1 - Jaccard similarity of top GENIE3 edges",
      x = "MDS1",
      y = "MDS2"
    )
  
  ggsave(
    filename = file.path(OUT_DIR, paste0("MDS_", tissue, "_species_distance.png")),
    plot = p, width = 6.5, height = 5.5, dpi = 300
  )
}

# create an MDS multiplot for manuscript
library(png)
library(grid)
library(gridExtra)

# Read the PNG images
img_B <- readPNG("/Users/tarangmehta/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/TGAC/Technology_Development/ATAC-Seq/ATAC_manuscript/Files/genie3/comparative_GRN/MDS_B_species_distance.png")
img_E <- readPNG("/Users/tarangmehta/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/TGAC/Technology_Development/ATAC-Seq/ATAC_manuscript/Files/genie3/comparative_GRN/MDS_E_species_distance.png")
img_L <- readPNG("/Users/tarangmehta/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/TGAC/Technology_Development/ATAC-Seq/ATAC_manuscript/Files/genie3/comparative_GRN/MDS_L_species_distance.png")
img_T <- readPNG("/Users/tarangmehta/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/TGAC/Technology_Development/ATAC-Seq/ATAC_manuscript/Files/genie3/comparative_GRN/MDS_T_species_distance.png")

# Convert to raster grobs (no interpolation to preserve plot sharpness)
g1 <- rasterGrob(img_B, interpolate = FALSE)
g2 <- rasterGrob(img_E, interpolate = FALSE)
g3 <- rasterGrob(img_L, interpolate = FALSE)
g4 <- rasterGrob(img_T, interpolate = FALSE)

# Save 2x2 multiplot as PNG (adjust width/height/res as needed)
png("/Users/tarangmehta/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/TGAC/Technology_Development/ATAC-Seq/ATAC_manuscript/Files/genie3/comparative_GRN/MDS_2x2_multiplot.png", width = 12, height = 10, units = "in", res = 300)
grid.arrange(g1, g2, g3, g4, ncol = 2, nrow = 2)
dev.off()

# -----------------------------
# Read orthologs
# -----------------------------
cat("Reading orthologs...\n")
orthologs <- fread(
  ORTHOLOGS_FILE,
  header = TRUE,
  fill = TRUE,
  na.strings = c("", "NA", "NULL")
)

orthologs <- orthologs %>%
  mutate(
    Orthogroup = as.character(Orthogroup),
    Unified_GeneSymbol = as.character(Unified_GeneSymbol)
  ) %>%
  distinct(Orthogroup, .keep_all = TRUE) %>%
  as.data.table()

# -----------------------------
# Rewired edges per tissue
# -----------------------------
rewired_all <- list()
tissue_rewired_plots <- list()

for (tissue in tissue_codes) {
  tissue_nets <- network_list[purrr::map_chr(network_list, "tissue") == tissue]
  if (length(tissue_nets) < 2) next
  
  species_names <- purrr::map_chr(tissue_nets, "species")
  
  combined <- purrr::map_dfr(tissue_nets, function(net) {
    net$top %>%
      select(Regulator, Target, Weight, Rank, Edge) %>%
      mutate(Species = net$species, Tissue = net$tissue)
  })
  
  edge_presence <- combined %>%
    group_by(Tissue, Edge, Regulator, Target) %>%
    summarise(
      n_species_present = n_distinct(Species),
      species_present = paste(sort(unique(Species)), collapse = ","),
      mean_weight = mean(Weight, na.rm = TRUE),
      .groups = "drop"
    )
  
  write.table(
    edge_presence,
    file.path(OUT_DIR, paste0("Edge_presence_summary_", tissue, ".tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  for (sp in species_names) {
    sp_edges <- combined %>% filter(Species == sp)
    
    others <- combined %>%
      filter(Species != sp) %>%
      group_by(Edge) %>%
      summarise(
        other_species_count = n_distinct(Species),
        other_mean_weight = mean(Weight, na.rm = TRUE),
        other_best_rank = min(Rank, na.rm = TRUE),
        .groups = "drop"
      )
    
    rewired <- sp_edges %>%
      left_join(edge_presence, by = c("Tissue", "Edge", "Regulator", "Target")) %>%
      left_join(others, by = "Edge") %>%
      mutate(
        other_species_count = ifelse(is.na(other_species_count), 0, other_species_count),
        other_mean_weight   = ifelse(is.na(other_mean_weight), 0, other_mean_weight),
        other_best_rank     = ifelse(is.na(other_best_rank), Inf, other_best_rank),
        unique_to_species   = n_species_present == 1,
        weight_delta        = Weight - other_mean_weight,
        rank_advantage      = ifelse(is.finite(other_best_rank), other_best_rank - Rank, Rank),
        rewiring_score      = ifelse(unique_to_species, 5, 0) +
          safe_scale(weight_delta) +
          safe_scale(rank_advantage)
      ) %>%
      arrange(desc(rewiring_score))
    
    # Orthogroup -> symbol mapping
    reg_map <- map_oma_to_symbol(rewired$Regulator, orthologs)
    tgt_map <- map_oma_to_symbol(rewired$Target, orthologs)
    
    rewired <- rewired %>%
      mutate(
        Regulator_Orthogroup = Regulator,
        Target_Orthogroup = Target,
        Regulator_UnifiedGS = reg_map$Unified_GeneSymbol,
        Target_UnifiedGS = tgt_map$Unified_GeneSymbol,
        Plot_Regulator = ifelse(!is.na(Regulator_UnifiedGS) & Regulator_UnifiedGS != "",
                                Regulator_UnifiedGS, Regulator),
        Plot_Target = ifelse(!is.na(Target_UnifiedGS) & Target_UnifiedGS != "",
                             Target_UnifiedGS, Target),
        Plot_Label = paste(Plot_Regulator, "\u2192", Plot_Target)
      )
    
    rewired_top <- head(rewired, min(TOP_REWIRED, nrow(rewired)))
    rewired_all[[paste(tissue, sp, sep = "_")]] <- rewired_top
    
    write.table(
      rewired,
      file.path(OUT_DIR, paste0(sp, "_", tissue, "_rewiring_all_edges.tsv")),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    
    write.table(
      rewired_top,
      file.path(OUT_DIR, paste0(sp, "_", tissue, "_top", TOP_REWIRED, "_rewired_edges.tsv")),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    
    plot_df <- rewired_top %>%
      arrange(rewiring_score) %>%
      mutate(Plot_Label = factor(Plot_Label, levels = Plot_Label))
    
    p <- ggplot(plot_df, aes(x = rewiring_score, y = Plot_Label, fill = unique_to_species)) +
      geom_col() +
      scale_fill_manual(values = c("FALSE" = "#80b1d3", "TRUE" = "#e41a1c")) +
      theme_bw(base_size = 11) +
      labs(
        title = paste("Top rewired edges:", sp, "-", tissue_names[[tissue]]),
        subtitle = "Labels use Unified_GeneSymbol when available",
        x = "Rewiring score",
        y = "TF \u2192 target",
        fill = "Unique edge"
      )
    
    ggsave(
      filename = file.path(OUT_DIR, paste0(sp, "_", tissue, "_top_rewired_edges.png")),
      plot = p, width = 8, height = 10, dpi = 300
    )
    
    tissue_rewired_plots[[paste(sp, tissue, sep = "_")]] <- p
  }
  
  # Rewiring summary with zero-count categories retained
  tissue_summary <- edge_presence %>%
    mutate(
      edge_class = case_when(
        n_species_present == 1 ~ "Species-specific",
        n_species_present >= MIN_SHARED_SPECIES & n_species_present < length(species_names) ~ "Shared in subset",
        n_species_present == length(species_names) ~ "Shared in all species",
        TRUE ~ "Other"
      )
    ) %>%
    count(edge_class) %>%
    filter(edge_class != "Other") %>%
    right_join(
      data.frame(edge_class = c("Species-specific", "Shared in subset", "Shared in all species")),
      by = "edge_class"
    ) %>%
    mutate(n = ifelse(is.na(n), 0, n)) %>%
    mutate(
      edge_class = factor(edge_class,
                          levels = c("Species-specific", "Shared in subset", "Shared in all species")),
      n_label = paste("n =", format(n, big.mark = ","))
    )
  
  write.table(
    tissue_summary,
    file.path(OUT_DIR, paste0("Rewiring_summary_", tissue, ".tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  p_sum <- ggplot(tissue_summary, aes(x = edge_class, y = n, fill = edge_class)) +
    geom_col(width = 0.75) +
    geom_text(aes(label = n_label),
              vjust = ifelse(tissue_summary$n > 0, -0.3, 0.5),
              size = 4) +
    scale_fill_manual(values = c(
      "Species-specific" = "#e63946",
      "Shared in subset" = "#f4a261",
      "Shared in all species" = "#2a9d8f"
    )) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none") +
    labs(
      title = tissue_names[[tissue]],
      x = NULL,
      y = "Number of top edges"
    )
  
  ggsave(
    filename = file.path(OUT_DIR, paste0("Rewiring_summary_", tissue, ".png")),
    plot = p_sum, width = 6.5, height = 5.5, dpi = 300
  )
}

## create stacked barplot of rewiring distribution for main figures - GENIE3 only

library(tidyverse)

base_dir <- "/Users/tarangmehta/Library/CloudStorage/OneDrive-TheUniversityOfLiverpool/TGAC/Technology_Development/ATAC-Seq/ATAC_manuscript/Files/genie3/comparative_GRN"

file_stems   <- c("B", "E", "L", "T")
tissue_names <- c("Brain", "Retina", "Liver", "Testis")
files <- file.path(base_dir, paste0("Rewiring_summary_", file_stems, ".tsv"))

label_threshold <- 0.0003   # 0.03%; change to 0.01 for 1%

plot_data <- purrr::map2_df(files, tissue_names, function(file, tissue) {
  d <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  d$tissue <- tissue
  d
}) |>
  mutate(
    n = as.numeric(n),
    tissue = factor(tissue, levels = c("Brain", "Retina", "Liver", "Testis")),
    edge_class = factor(
      edge_class,
      levels = c("Shared in all species", "Shared in subset", "Species-specific")
    )
  ) |>
  group_by(tissue) |>
  arrange(edge_class, .by_group = TRUE) |>
  mutate(
    total = sum(n),
    prop = n / total,
    pct_label = scales::percent(prop, accuracy = 0.1),
    n_label_clean = paste0("n = ", format(n, big.mark = ",")),
    label = ifelse(
      prop >= label_threshold,
      paste0(pct_label, "\n", n_label_clean),
      NA
    )
  ) |>
  ungroup()

p <- ggplot(plot_data, aes(x = tissue, y = prop, fill = edge_class)) +
  geom_col(color = "white", linewidth = 0.3, width = 0.75) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 4.5,
    lineheight = 0.95,
    color = "black",
    na.rm = TRUE
  ) +
  coord_flip() +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = c(
      "Shared in all species" = "#ffe6e6",
      "Shared in subset" = "#D55E00",
      "Species-specific" = "#009E73"
    )
  ) +
  labs(
    x = NULL,
    y = "Percentage of top 12,000 non-redundant rewired GENIE3 edges",
    fill = NULL
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_blank(),
    axis.text.y = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_blank(),
    legend.text = element_text(size = 16),
    legend.position = "bottom",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

# Add these theme tweaks FIRST (before ggsave) - ensures text beyond margins are visible
p_final <- p +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.05)),  # Add 5% right padding
    breaks = c(0, 0.25, 0.5, 0.75, 1.0)    # Force 100% to show
  ) +
  theme(
    plot.margin = margin(15, 30, 15, 15),  # Extra right margin
    text = element_text(size = 20)         # Bump up all text
  )

# Export settings for publication
ggsave(
  filename = file.path("/Users/tarangmehta/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/TGAC/Technology_Development/ATAC-Seq/ATAC_manuscript/Figures/fig6_or_supp/6a_genie3_rewiringcategories.png"),
  plot = p_final,
  width = 26,      # 18cm = ~7 inches (perfect for 4 bars)
  height = 10,      # 8cm = ~3.15 inches (fits 4 bars + legend)
  units = "cm",    # Explicitly specify cm
  dpi = 300,
  bg = "white"
)

# -----------------------------
# Master summary table
# -----------------------------
summary_df <- file_table %>%
  mutate(Network_ID = paste(Species, Tissue, sep = "_")) %>%
  rowwise() %>%
  mutate(
    Total_edges = nrow(network_list[[Network_ID]]$full),
    Top_edges_used = nrow(network_list[[Network_ID]]$top),
    Unique_regulators = dplyr::n_distinct(network_list[[Network_ID]]$top$Regulator),
    Unique_targets = dplyr::n_distinct(network_list[[Network_ID]]$top$Target)
  ) %>%
  ungroup()

write.table(
  summary_df,
  file.path(OUT_DIR, "comparative_GRN_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# -----------------------------
# Jaccard multiplot
# -----------------------------
jaccard_plot_list <- list()

for (tissue in tissue_codes) {
  f <- file.path(OUT_DIR, paste0("Jaccard_", tissue, "_heatmap.png"))
  if (file.exists(f)) {
    img <- png::readPNG(f)
    jaccard_plot_list[[tissue_names[[tissue]]]] <- ggplot() +
      annotation_raster(img, -Inf, Inf, -Inf, Inf) +
      theme_void() +
      labs(title = tissue_names[[tissue]]) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
      )
  }
}

if (length(jaccard_plot_list) > 0) {
  multi_jaccard <- (jaccard_plot_list[["Brain"]] | jaccard_plot_list[["Retina"]]) /
    (jaccard_plot_list[["Liver"]] | jaccard_plot_list[["Testis"]]) +
    plot_annotation(
      title = "Tissue-specific Jaccard similarity of top 12,000 GENIE3 edges across species",
      subtitle = paste("Top", format(TOP_N_EDGES, big.mark = ","), "edges per network")
    )
  
  ggsave(
    filename = file.path(OUT_DIR, "Jaccard_all_tissues_labeled.png"),
    plot = multi_jaccard, width = 14, height = 12, dpi = 300, bg = "white"
  )
}

# -----------------------------
# Rewired-edge multiplots per tissue
# -----------------------------
for (tissue in tissue_codes) {
  keys <- names(tissue_rewired_plots)[str_detect(names(tissue_rewired_plots), paste0("_", tissue, "$"))]
  tissue_plots <- tissue_rewired_plots[keys]
  
  if (length(tissue_plots) >= 2) {
    multi_rewired <- wrap_plots(tissue_plots, ncol = 2) +
      plot_annotation(
        title = paste("Top rewired GENIE3 edges -", tissue_names[[tissue]]),
        subtitle = "Unified gene symbols used in plot labels where available"
      )
    
    ggsave(
      filename = file.path(OUT_DIR, paste0("Rewired_edges_", tissue, "_multiplot.png")),
      plot = multi_rewired, width = 16, height = 30, dpi = 300, limitsize = FALSE
    )
  }
}

# -----------------------------
# Rewiring summary multiplot
# -----------------------------
summary_plot_list <- list()

for (tissue in tissue_codes) {
  f <- file.path(OUT_DIR, paste0("Rewiring_summary_", tissue, ".png"))
  if (file.exists(f)) {
    img <- png::readPNG(f)
    summary_plot_list[[tissue_names[[tissue]]]] <- ggplot() +
      annotation_raster(img, -Inf, Inf, -Inf, Inf) +
      theme_void() +
      labs(title = tissue_names[[tissue]]) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
      )
  }
}

if (length(summary_plot_list) > 0) {
  multi_summary <- (summary_plot_list[["Brain"]] | summary_plot_list[["Retina"]]) /
    (summary_plot_list[["Liver"]] | summary_plot_list[["Testis"]]) +
    plot_annotation(
      title = "GENIE3 edge conservation vs rewiring across tissues",
      subtitle = "Includes zero-count categories and n labels on bars"
    )
  
  ggsave(
    filename = file.path(OUT_DIR, "Rewiring_summary_all_tissues_labeled.png"),
    plot = multi_summary, width = 14, height = 12, dpi = 300, bg = "white"
  )
}

cat("\nDone.\nOutputs written to:\n", OUT_DIR, "\n")
cat("\nKey outputs:\n")
cat("- Jaccard_[B/E/L/T]_heatmap.png\n")
cat("- Jaccard_all_tissues_labeled.png\n")
cat("- *_rewiring_all_edges.tsv\n")
cat("- *_top100_rewired_edges.tsv\n")
cat("- *_top_rewired_edges.png\n")
cat("- Rewired_edges_[B/E/L/T]_multiplot.png\n")
cat("- Rewiring_summary_[B/E/L/T].png\n")
cat("- Rewiring_summary_all_tissues_labeled.png\n")
cat("- comparative_GRN_summary.tsv\n")

#############################################################################################################

## to do next, filter by TF footprints of Regulator and Target, then replot and compare - what changes
# Mostly report the new results but comment on whether proportion of '12k top edges' in Rewiring_summary_*.png and "Jaccard_all_tissues_labeled.png" plots changes from expression only to expression + footprints
# Is there added inference by including TF footprints (likely)

# ================================================
# Complete GENIE3 Comparative GRN Analysis - NOW USING COMPOSITE WEIGHT BASED ON BOTH GENE EXPRESSION AND TF FOOTPRINTS (CompositeWeight=GENIE3_Weight×MotifBitScorenormalized)
# - Reads *_GENIE3_links.composite.tsv
# - Keeps the top n edges per network ranked by CompositeWeight
# - Computes tissue-specific Jaccard matrices
# - Plots heatmaps
# - Extracts top rewired edges using CompositeWeight-based measures
# - Maps Orthogroup IDs to Unified_GeneSymbol
# - Builds publication-style multiplots
# ================================================

# Composite Weight explained:
# GENIE3 weight = predicted regulatory strength
# Motif bit-score = binding site quality
# Product = regulatory strength × binding evidence (both must be strong)
# 
# Normalization first - Motif bit-scores need scaling to match GENIE3's 0-1 range so this uses:
# MotifBitScorenorm=BitScore−BitScoremin/⁡BitScoremax⁡−BitScoremin⁡
# All above is done in: /Users/tarangmehta/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/github/ATAC_bioinformatics/ATAC_Bioinf_ManuscriptAnalysis.sh (see LINE 8992 ONWARDS)


# -----------------------------
# Packages
# -----------------------------
required_packages <- c(
  "data.table", "dplyr", "tidyr", "ggplot2", "patchwork",
  "pheatmap", "stringr", "purrr", "readr", "png", "tibble"
)

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

# -----------------------------
# Paths and parameters
# -----------------------------
DATA_DIR <- "/Users/tarangmehta/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/TGAC/Technology_Development/ATAC-Seq/ATAC_manuscript/Files/genie3/genie3_composite"
OUT_DIR  <- file.path(DATA_DIR, "comparative_GRN")
ORTHOLOGS_FILE <- "/Users/tarangmehta/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/TGAC/Technology_Development/ATAC-Seq/ATAC_manuscript/Files/OrthologousGroups_ENS_GeneSymbols_6cich.txt"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Compare only same tissues across species
tissue_codes <- c("B", "E", "L", "T")

# Tissue names
tissue_names <- c(
  "B" = "Brain",
  "E" = "Retina",
  "L" = "Liver",
  "T" = "Testis"
)

# Keep top N edges from each network
TOP_N_EDGES <- 12000

# Number of top rewired edges to report
TOP_REWIRED <- 100

# Optional minimum number of species an edge must appear in to be called shared
MIN_SHARED_SPECIES <- 2

# -----------------------------
# Helper functions
# -----------------------------
get_edge_id <- function(df) {
  paste(df$Regulator, df$Target, sep = "->")
}

jaccard_index <- function(edges1, edges2) {
  inter <- length(intersect(edges1, edges2))
  union <- length(union(edges1, edges2))
  if (union == 0) return(NA_real_)
  inter / union
}

safe_read_links <- function(file) {
  df <- fread(file, header = TRUE, data.table = FALSE)
  
  expected_cols <- c(
    "Regulator", "Target", "Weight",
    "Motif_BitScore", "MotifBitScore_norm", "CompositeWeight"
  )
  
  if (ncol(df) < 6) {
    stop(paste("File has fewer than 6 columns:", file))
  }
  
  colnames(df)[1:6] <- expected_cols
  
  df <- df %>%
    filter(
      !is.na(Regulator),
      !is.na(Target),
      !is.na(Weight),
      !is.na(Motif_BitScore),
      !is.na(MotifBitScore_norm),
      !is.na(CompositeWeight)
    ) %>%
    mutate(
      Regulator = as.character(Regulator),
      Target = as.character(Target),
      Weight = as.numeric(Weight),
      Motif_BitScore = as.numeric(Motif_BitScore),
      MotifBitScore_norm = as.numeric(MotifBitScore_norm),
      CompositeWeight = as.numeric(CompositeWeight)
    ) %>%
    filter(
      !is.na(Weight),
      !is.na(Motif_BitScore),
      !is.na(MotifBitScore_norm),
      !is.na(CompositeWeight)
    ) %>%
    arrange(desc(CompositeWeight), desc(Weight)) %>%
    mutate(
      Rank = row_number(),
      Edge = paste(Regulator, Target, sep = "->")
    )
  
  df
}

map_oma_to_symbol <- function(gene_col, orthologs_dt) {
  out <- data.table(
    Original_Gene = as.character(gene_col),
    Unified_GeneSymbol = NA_character_,
    Mapping_Success = FALSE
  )
  
  idx <- match(out$Original_Gene, orthologs_dt$Orthogroup)
  hit <- !is.na(idx)
  
  out$Unified_GeneSymbol[hit] <- orthologs_dt$Unified_GeneSymbol[idx[hit]]
  out$Mapping_Success[hit] <- !is.na(out$Unified_GeneSymbol[hit]) & out$Unified_GeneSymbol[hit] != ""
  
  out
}

safe_scale <- function(x) {
  if (all(is.na(x))) return(rep(0, length(x)))
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  as.numeric(scale(x))
}

# -----------------------------
# Find all GENIE3 composite link files
# -----------------------------
all_files <- list.files(
  DATA_DIR,
  pattern = "^[A-Za-z]+_[BELT]_GENIE3_links.bitscore\\.composite\\.tsv$",
  full.names = TRUE
)

if (length(all_files) == 0) {
  stop("No *_GENIE3_links.composite.tsv files found in DATA_DIR")
}

file_table <- tibble(
  file = all_files,
  base = basename(all_files)
) %>%
  mutate(
    Species = str_extract(base, "^[A-Za-z]+"),
    Tissue  = str_extract(base, "(?<=_)[BELT](?=_GENIE3_links.bitscore\\.composite\\.tsv$)")
  ) %>%
  arrange(Tissue, Species)

write.table(
  file_table,
  file.path(OUT_DIR, "detected_GENIE3_composite_link_files.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# -----------------------------
# Read all networks
# -----------------------------
network_list <- list()

for (i in seq_len(nrow(file_table))) {
  f <- file_table$file[i]
  sp <- file_table$Species[i]
  ts <- file_table$Tissue[i]
  
  cat("Reading:", basename(f), "\n")
  df <- safe_read_links(f)
  df_top <- head(df, min(TOP_N_EDGES, nrow(df)))
  
  network_list[[paste(sp, ts, sep = "_")]] <- list(
    species = sp,
    tissue = ts,
    full = df,      # ALL edges
    top = df_top    # Top 12k
  )
}

# -----------------------------
# Tissue-specific Jaccard matrices
# -----------------------------
jaccard_summary <- list()

for (tissue in tissue_codes) {
  tissue_nets <- network_list[purrr::map_chr(network_list, "tissue") == tissue]
  
  if (length(tissue_nets) < 2) {
    cat("Skipping tissue", tissue, "- fewer than 2 species available\n")
    next
  }
  
  species_names <- purrr::map_chr(tissue_nets, "species")
  edge_sets <- purrr::map(tissue_nets, ~ .x$top$Edge)
  
  mat <- matrix(
    NA_real_,
    nrow = length(species_names),
    ncol = length(species_names),
    dimnames = list(species_names, species_names)
  )
  
  for (i in seq_along(species_names)) {
    for (j in seq_along(species_names)) {
      mat[i, j] <- jaccard_index(edge_sets[[i]], edge_sets[[j]])
    }
  }
  
  jaccard_summary[[tissue]] <- as.data.frame(mat)
  
  write.table(
    mat,
    file.path(OUT_DIR, paste0("Jaccard_", tissue, "_species_matrix.tsv")),
    sep = "\t", quote = FALSE, col.names = NA
  )
  
  pheatmap::pheatmap(
    mat,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    display_numbers = TRUE,
    number_format = "%.2f",
    color = colorRampPalette(c("white", "#fee08b", "#d73027"))(100),
    border_color = NA,
    main = paste("Jaccard similarity of top", TOP_N_EDGES, "composite-ranked GENIE3 edges -", tissue_names[[tissue]]),
    filename = file.path(OUT_DIR, paste0("Jaccard_", tissue, "_heatmap.png")),
    width = 7,
    height = 6
  )
  
  long_df <- as.data.frame(as.table(mat)) %>%
    rename(Species1 = Var1, Species2 = Var2, Jaccard = Freq) %>%
    arrange(desc(Jaccard))
  
  write.table(
    long_df,
    file.path(OUT_DIR, paste0("Jaccard_", tissue, "_long.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

# -----------------------------
# Distance / MDS plots per tissue
# -----------------------------
for (tissue in names(jaccard_summary)) {
  mat <- as.matrix(jaccard_summary[[tissue]])
  dist_mat <- 1 - mat
  diag(dist_mat) <- 0
  
  mds <- cmdscale(as.dist(dist_mat), k = 2, eig = TRUE)
  mds_df <- data.frame(
    Species = rownames(mat),
    Dim1 = mds$points[, 1],
    Dim2 = mds$points[, 2]
  )
  
  write.table(
    mds_df,
    file.path(OUT_DIR, paste0("MDS_", tissue, "_species_positions.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  p <- ggplot(mds_df, aes(Dim1, Dim2, label = Species)) +
    geom_point(size = 4, color = "#2c7fb8") +
    geom_text(vjust = -0.8, size = 4) +
    theme_bw(base_size = 12) +
    labs(
      title = paste("Species distance in GRN space -", tissue_names[[tissue]]),
      subtitle = "Distance = 1 - Jaccard similarity of top composite-ranked GENIE3 edges",
      x = "MDS1",
      y = "MDS2"
    )
  
  ggsave(
    filename = file.path(OUT_DIR, paste0("MDS_", tissue, "_species_distance.png")),
    plot = p, width = 6.5, height = 5.5, dpi = 300
  )
}

# create an MDS multiplot for manuscript
library(png)
library(grid)
library(gridExtra)

# Read the PNG images
img_B2 <- readPNG("/Users/tarangmehta/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/TGAC/Technology_Development/ATAC-Seq/ATAC_manuscript/Files/genie3/genie3_composite/comparative_GRN/MDS_B_species_distance.png")
img_E2 <- readPNG("/Users/tarangmehta/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/TGAC/Technology_Development/ATAC-Seq/ATAC_manuscript/Files/genie3/genie3_composite/comparative_GRN/MDS_E_species_distance.png")
img_L2 <- readPNG("/Users/tarangmehta/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/TGAC/Technology_Development/ATAC-Seq/ATAC_manuscript/Files/genie3/genie3_composite/comparative_GRN/MDS_L_species_distance.png")
img_T2 <- readPNG("/Users/tarangmehta/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/TGAC/Technology_Development/ATAC-Seq/ATAC_manuscript/Files/genie3/genie3_composite/comparative_GRN/MDS_T_species_distance.png")

# Convert to raster grobs (no interpolation to preserve plot sharpness)
g1a <- rasterGrob(img_B2, interpolate = FALSE)
g2a <- rasterGrob(img_E2, interpolate = FALSE)
g3a <- rasterGrob(img_L2, interpolate = FALSE)
g4a <- rasterGrob(img_T2, interpolate = FALSE)

# Save 2x2 multiplot as PNG (adjust width/height/res as needed)
png("/Users/tarangmehta/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/TGAC/Technology_Development/ATAC-Seq/ATAC_manuscript/Files/genie3/genie3_composite/comparative_GRN/MDS_2x2_multiplot.png", width = 12, height = 10, units = "in", res = 300)
grid.arrange(g1a, g2a, g3a, g4a, ncol = 2, nrow = 2)
dev.off()

# -----------------------------
# Read orthologs
# -----------------------------
cat("Reading orthologs...\n")
orthologs <- fread(
  ORTHOLOGS_FILE,
  header = TRUE,
  fill = TRUE,
  na.strings = c("", "NA", "NULL")
)

orthologs <- orthologs %>%
  mutate(
    Orthogroup = as.character(Orthogroup),
    Unified_GeneSymbol = as.character(Unified_GeneSymbol)
  ) %>%
  distinct(Orthogroup, .keep_all = TRUE) %>%
  as.data.table()

# -----------------------------
# Rewired edges per tissue
# -----------------------------
rewired_all <- list()                          # top 100 rewired subset per species-tissue  
rewired_background_12k <- list()               # top 12k composite background per species-tissue
all_edges_rewired <- list()                    # ALL rewired edges per species-tissue (NEW)
tissue_rewired_plots <- list()

for (tissue in tissue_codes) {
  tissue_nets <- network_list[purrr::map_chr(network_list, "tissue") == tissue]
  if (length(tissue_nets) < 2) next
  
  species_names <- purrr::map_chr(tissue_nets, "species")
  
  combined <- purrr::map_dfr(tissue_nets, function(net) {
    net$top %>%  # Still use top 12k for combined (conservative)
      select(
        Regulator, Target, Weight, Motif_BitScore,
        MotifBitScore_norm, CompositeWeight, Rank, Edge
      ) %>%
      mutate(Species = net$species, Tissue = net$tissue)
  })
  
  edge_presence <- combined %>%
    group_by(Tissue, Edge, Regulator, Target) %>%
    summarise(
      n_species_present = n_distinct(Species),
      species_present = paste(sort(unique(Species)), collapse = ","),
      mean_weight = mean(Weight, na.rm = TRUE),
      mean_composite_weight = mean(CompositeWeight, na.rm = TRUE),
      mean_motif_bitscore = mean(Motif_BitScore, na.rm = TRUE),
      mean_motif_bitscore_norm = mean(MotifBitScore_norm, na.rm = TRUE),
      .groups = "drop"
    )
  
  write.table(
    edge_presence,
    file.path(OUT_DIR, paste0("Edge_presence_summary_", tissue, ".tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  for (sp in species_names) {
    sp_edges_top12k <- combined %>% filter(Species == sp)  # Top 12k background
    
    # ALL edges with Species/Tissue added (FIX)
    sp_edges_full <- network_list[[paste(sp, tissue, sep="_")]]$full %>%
      mutate(Species = sp, Tissue = tissue)  # ADD these columns
    
    # Save top 12k background (unchanged)
    rewired_background_12k[[paste(tissue, sp, sep = "_")]] <- sp_edges_top12k
    
    others <- combined %>%
      filter(Species != sp) %>%
      group_by(Edge) %>%
      summarise(
        other_species_count = n_distinct(Species),
        other_mean_weight = mean(Weight, na.rm = TRUE),
        other_mean_composite_weight = mean(CompositeWeight, na.rm = TRUE),
        other_mean_motif_bitscore = mean(Motif_BitScore, na.rm = TRUE),
        other_mean_motif_bitscore_norm = mean(MotifBitScore_norm, na.rm = TRUE),
        other_best_rank = min(Rank, na.rm = TRUE),
        .groups = "drop"
      )
    
    # Rewiring on FULL edges WITH Tissue/Species
    rewired_full <- sp_edges_full %>%
      left_join(edge_presence, by = c("Tissue", "Edge", "Regulator", "Target")) %>%
      left_join(others, by = "Edge") %>%
      mutate(
        other_species_count = ifelse(is.na(other_species_count), 0, other_species_count),
        other_mean_weight = ifelse(is.na(other_mean_weight), 0, other_mean_weight),
        other_mean_composite_weight = ifelse(is.na(other_mean_composite_weight), 0, other_mean_composite_weight),
        other_mean_motif_bitscore = ifelse(is.na(other_mean_motif_bitscore), 0, other_mean_motif_bitscore),
        other_mean_motif_bitscore_norm = ifelse(is.na(other_mean_motif_bitscore_norm), 0, other_mean_motif_bitscore_norm),
        other_best_rank = ifelse(is.na(other_best_rank), Inf, other_best_rank),
        unique_to_species = ifelse(is.na(n_species_present), TRUE, n_species_present == 1),  # Handle NA
        composite_weight_delta = CompositeWeight - coalesce(other_mean_composite_weight, 0),
        rank_advantage = ifelse(is.finite(other_best_rank), other_best_rank - Rank, Rank),
        rewiring_score = ifelse(unique_to_species, 5, 0) +
          safe_scale(composite_weight_delta) +
          safe_scale(rank_advantage)
      ) %>%
      arrange(desc(rewiring_score), desc(CompositeWeight), desc(Weight))
    
    # Gene symbol mapping (on full)
    reg_map <- map_oma_to_symbol(rewired_full$Regulator, orthologs)
    tgt_map <- map_oma_to_symbol(rewired_full$Target, orthologs)
    
    rewired_full <- rewired_full %>%
      mutate(
        Regulator_Orthogroup = Regulator,
        Target_Orthogroup = Target,
        Regulator_UnifiedGS = reg_map$Unified_GeneSymbol,
        Target_UnifiedGS = tgt_map$Unified_GeneSymbol,
        Plot_Regulator = ifelse(!is.na(Regulator_UnifiedGS) & Regulator_UnifiedGS != "", Regulator_UnifiedGS, Regulator),
        Plot_Target = ifelse(!is.na(Target_UnifiedGS) & Target_UnifiedGS != "", Target_UnifiedGS, Target),
        Plot_Label = paste(Plot_Regulator, "→", Plot_Target)
      )
    
    # Save ALL rewired edges as "_rewiring_all_edges.tsv" 
    all_edges_rewired[[paste(tissue, sp, sep = "_")]] <- rewired_full
    write.table(
      rewired_full,
      file.path(OUT_DIR, paste0(sp, "_", tissue, "_rewiring_all_edges.tsv")),  # Now ALL edges
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    
    # Top 100 rewired subset (unchanged)
    rewired_top <- head(rewired_full, min(TOP_REWIRED, nrow(rewired_full)))
    rewired_all[[paste(tissue, sp, sep = "_")]] <- rewired_top
    
    write.table(
      sp_edges_top12k,  # Top 12k unchanged
      file.path(OUT_DIR, paste0(sp, "_", tissue, "_background_top12k_edges.tsv")),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    
    write.table(
      rewired_top,
      file.path(OUT_DIR, paste0(sp, "_", tissue, "_top", TOP_REWIRED, "_rewired_edges.tsv")),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    
    # Plot top 100 (unchanged)
    plot_df <- rewired_top %>%
      arrange(rewiring_score) %>%
      mutate(
        Plot_Label_raw = paste(Plot_Regulator, "→", Plot_Target),
        Plot_Label = make.unique(Plot_Label_raw),
        Plot_Label = factor(Plot_Label, levels = rev(Plot_Label))
      )
    
    p <- ggplot(plot_df, aes(x = rewiring_score, y = Plot_Label, fill = unique_to_species)) +
      geom_col() +
      scale_fill_manual(values = c("FALSE" = "#80b1d3", "TRUE" = "#e41a1c")) +
      theme_bw(base_size = 11) +
      labs(
        title = paste("Top rewired edges:", sp, "-", tissue_names[[tissue]]),
        subtitle = "CompositeWeight-based rewiring; labels use Unified_GeneSymbol when available",
        x = "Rewiring score",
        y = "TF → target",
        fill = "Unique edge"
      )
    
    ggsave(
      filename = file.path(OUT_DIR, paste0(sp, "_", tissue, "_top_rewired_edges.png")),
      plot = p, width = 12, height = 14, dpi = 300
    )
    
    tissue_rewired_plots[[paste(sp, tissue, sep = "_")]] <- p
  }
  
  # Tissue summary/plot (unchanged)
  tissue_summary <- edge_presence %>%
    mutate(
      edge_class = case_when(
        n_species_present == 1 ~ "Species-specific",
        n_species_present >= MIN_SHARED_SPECIES & n_species_present < length(species_names) ~ "Shared in subset",
        n_species_present == length(species_names) ~ "Shared in all species",
        TRUE ~ "Other"
      )
    ) %>%
    count(edge_class) %>%
    filter(edge_class != "Other") %>%
    right_join(
      data.frame(edge_class = c("Species-specific", "Shared in subset", "Shared in all species")),
      by = "edge_class"
    ) %>%
    mutate(n = ifelse(is.na(n), 0, n)) %>%
    mutate(
      edge_class = factor(
        edge_class,
        levels = c("Species-specific", "Shared in subset", "Shared in all species")
      ),
      n_label = paste("n =", format(n, big.mark = ","))
    )
  
  write.table(
    tissue_summary,
    file.path(OUT_DIR, paste0("Rewiring_summary_", tissue, ".tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  p_sum <- ggplot(tissue_summary, aes(x = edge_class, y = n, fill = edge_class)) +
    geom_col(width = 0.75) +
    geom_text(
      aes(label = n_label),
      vjust = ifelse(tissue_summary$n > 0, -0.3, 0.5),
      size = 4
    ) +
    scale_fill_manual(values = c(
      "Species-specific" = "#e63946",
      "Shared in subset" = "#f4a261",
      "Shared in all species" = "#2a9d8f"
    )) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none") +
    labs(
      title = tissue_names[[tissue]],
      x = NULL,
      y = "Number of top edges"
    )
  
  ggsave(
    filename = file.path(OUT_DIR, paste0("Rewiring_summary_", tissue, ".png")),
    plot = p_sum, width = 6.5, height = 5.5, dpi = 300
  )
}

# # Save lists for downstream
# saveRDS(all_edges_rewired, file.path(OUT_DIR, "all_edges_rewired.rds"))
# saveRDS(rewired_all, file.path(OUT_DIR, "rewired_top100.rds"))
# saveRDS(rewired_background_12k, file.path(OUT_DIR, "rewired_background_top12k.rds"))

## create stacked barplot of rewiring distribution for main figures - GENIE3 only

library(tidyverse)

base_dir <- "/Users/tarangmehta/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/TGAC/Technology_Development/ATAC-Seq/ATAC_manuscript/Files/genie3/genie3_composite/comparative_GRN"

file_stems   <- c("B", "E", "L", "T")
tissue_names <- c("Brain", "Retina", "Liver", "Testis")
files <- file.path(base_dir, paste0("Rewiring_summary_", file_stems, ".tsv"))

label_threshold <- 0.01   # 0.01 for 1%

plot_data <- purrr::map2_df(files, tissue_names, function(file, tissue) {
  d <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  d$tissue <- tissue
  d
}) |>
  mutate(
    n = as.numeric(n),
    tissue = factor(tissue, levels = c("Brain", "Retina", "Liver", "Testis")),
    edge_class = factor(
      edge_class,
      levels = c("Shared in all species", "Shared in subset", "Species-specific")
    )
  ) |>
  group_by(tissue) |>
  arrange(edge_class, .by_group = TRUE) |>
  mutate(
    total = sum(n),
    prop = n / total,
    pct_label = scales::percent(prop, accuracy = 0.1),
    n_label_clean = paste0("n = ", format(n, big.mark = ",")),
    label = ifelse(
      prop >= label_threshold,
      paste0(pct_label, "\n", n_label_clean),
      NA
    )
  ) |>
  ungroup()

p <- ggplot(plot_data, aes(x = tissue, y = prop, fill = edge_class)) +
  geom_col(color = "white", linewidth = 0.3, width = 0.75) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 4.5,
    lineheight = 0.95,
    color = "black",
    na.rm = TRUE
  ) +
  coord_flip() +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = c(
      "Shared in all species" = "#ffe6e6",
      "Shared in subset" = "#D55E00",
      "Species-specific" = "#009E73"
    )
  ) +
  labs(
    x = NULL,
    y = "Percentage of top 12,000 non-redundant rewired motif-supported GENIE3 edges",
    fill = NULL
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_blank(),
    axis.text.y = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_blank(),
    legend.text = element_text(size = 16),
    legend.position = "bottom",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

# Add these theme tweaks FIRST (before ggsave) - ensures text beyond margins are visible
p_final <- p +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.05)),  # Add 5% right padding
    breaks = c(0, 0.25, 0.5, 0.75, 1.0)    # Force 100% to show
  ) +
  theme(
    plot.margin = margin(15, 30, 15, 15),  # Extra right margin
    text = element_text(size = 20)         # Bump up all text
  )

# Export settings for publication
ggsave(
  filename = file.path("/Users/tarangmehta/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/TGAC/Technology_Development/ATAC-Seq/ATAC_manuscript/Figures/fig6_or_supp/6b_motifGENIE3_rewiringcategories.png"),
  plot = p_final,
  width = 26,      # 18cm = ~7 inches (perfect for 4 bars)
  height = 10,      # 8cm = ~3.15 inches (fits 4 bars + legend)
  units = "cm",    # Explicitly specify cm
  dpi = 300,
  bg = "white"
)

# -----------------------------
# Master summary table
# -----------------------------
summary_df <- file_table %>%
  mutate(Network_ID = paste(Species, Tissue, sep = "_")) %>%
  rowwise() %>%
  mutate(
    Total_edges = nrow(network_list[[Network_ID]]$full),
    Top_edges_used = nrow(network_list[[Network_ID]]$top),
    Unique_regulators = dplyr::n_distinct(network_list[[Network_ID]]$top$Regulator),
    Unique_targets = dplyr::n_distinct(network_list[[Network_ID]]$top$Target)
  ) %>%
  ungroup()

write.table(
  summary_df,
  file.path(OUT_DIR, "comparative_GRN_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# -----------------------------
# Jaccard multiplot
# -----------------------------
jaccard_plot_list <- list()

for (tissue in tissue_codes) {
  f <- file.path(OUT_DIR, paste0("Jaccard_", tissue, "_heatmap.png"))
  if (file.exists(f)) {
    img <- png::readPNG(f)
    jaccard_plot_list[[tissue_names[[tissue]]]] <- ggplot() +
      annotation_raster(img, -Inf, Inf, -Inf, Inf) +
      theme_void() +
      labs(title = tissue_names[[tissue]]) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
      )
  }
}

if (length(jaccard_plot_list) > 0 &&
    all(c("Brain", "Retina", "Liver", "Testis") %in% names(jaccard_plot_list))) {
  multi_jaccard <- (jaccard_plot_list[["Brain"]] | jaccard_plot_list[["Retina"]]) /
    (jaccard_plot_list[["Liver"]] | jaccard_plot_list[["Testis"]]) +
    plot_annotation(
      title = "Tissue-specific Jaccard similarity of top composite-ranked GENIE3 edges across species",
      subtitle = paste("Top", format(TOP_N_EDGES, big.mark = ","), "edges per network ranked by CompositeWeight")
    )
  
  ggsave(
    filename = file.path(OUT_DIR, "Jaccard_all_tissues_labeled.png"),
    plot = multi_jaccard, width = 14, height = 12, dpi = 300, bg = "white"
  )
}

# -----------------------------
# Rewired-edge multiplots per tissue
# -----------------------------
for (tissue in tissue_codes) {
  keys <- names(tissue_rewired_plots)[str_detect(names(tissue_rewired_plots), paste0("_", tissue, "$"))]
  tissue_plots <- tissue_rewired_plots[keys]
  
  if (length(tissue_plots) >= 2) {
    multi_rewired <- wrap_plots(tissue_plots, ncol = 2) +
      plot_annotation(
        title = paste("Top rewired GENIE3 edges -", tissue_names[[tissue]]),
        subtitle = "CompositeWeight-based rewiring; Unified gene symbols used in labels where available"
      )
    
    ggsave(
      filename = file.path(OUT_DIR, paste0("Rewired_edges_", tissue, "_multiplot.png")),
      plot = multi_rewired, width = 20, height = 30, dpi = 300, limitsize = FALSE
    )
  }
}

# -----------------------------
# Rewiring summary multiplot
# -----------------------------
summary_plot_list <- list()

for (tissue in tissue_codes) {
  f <- file.path(OUT_DIR, paste0("Rewiring_summary_", tissue, ".png"))
  if (file.exists(f)) {
    img <- png::readPNG(f)
    summary_plot_list[[tissue_names[[tissue]]]] <- ggplot() +
      annotation_raster(img, -Inf, Inf, -Inf, Inf) +
      theme_void() +
      labs(title = tissue_names[[tissue]]) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
      )
  }
}

if (length(summary_plot_list) > 0 &&
    all(c("Brain", "Retina", "Liver", "Testis") %in% names(summary_plot_list))) {
  multi_summary <- (summary_plot_list[["Brain"]] | summary_plot_list[["Retina"]]) /
    (summary_plot_list[["Liver"]] | summary_plot_list[["Testis"]]) +
    plot_annotation(
      title = "GENIE3 edge conservation vs rewiring across tissues",
      subtitle = "CompositeWeight-ranked edges; includes zero-count categories and n labels on bars"
    )
  
  ggsave(
    filename = file.path(OUT_DIR, "Rewiring_summary_all_tissues_labeled.png"),
    plot = multi_summary, width = 14, height = 12, dpi = 300, bg = "white"
  )
}

cat("\nDone.\nOutputs written to:\n", OUT_DIR, "\n")
cat("\nKey outputs:\n")
cat("- Jaccard_[B/E/L/T]_heatmap.png\n")
cat("- Jaccard_all_tissues_labeled.png\n")
cat("- *_rewiring_all_edges.tsv\n")
cat("- *_top100_rewired_edges.tsv\n")
cat("- *_top_rewired_edges.png\n")
cat("- Rewired_edges_[B/E/L/T]_multiplot.png\n")
cat("- Rewiring_summary_[B/E/L/T].png\n")
cat("- Rewiring_summary_all_tissues_labeled.png\n")
cat("- comparative_GRN_summary.tsv\n")

# -----------------------------
## HA-HE genes (high accessibility + high expression) are prime candidates for active regulatory elements.
# 3-part validation on ALL rewired edges (full *_rewiring_all_edges.tsv):
# - Enrichment: HA-HE categories overrepresented in top 12k rewired vs ALL rewired?
# - Score distributions: HA-HE edges higher composite/rewiring scores?
# - Top 50 HA-HE candidates by composite × rewiring
# -----------------------------

library(tidyverse)
library(ggdist)
library(patchwork)

OUT_DIR <- "/Users/tarangmehta/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/TGAC/Technology_Development/ATAC-Seq/ATAC_manuscript/Files/genie3/genie3_composite/comparative_GRN/AEanalysis"
AE_DIR <- OUT_DIR
GRN_DIR <- "/Users/tarangmehta/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/TGAC/Technology_Development/ATAC-Seq/ATAC_manuscript/Files/genie3/genie3_composite/comparative_GRN"
TOP_REWIRED <- 12000  # Top N for "foreground" enrichment test

dir.create(AE_DIR, showWarnings = FALSE, recursive = TRUE)

# 1. Load AE categories
ae_files <- list.files(AE_DIR, pattern = "_AEcateg.tsv$", full.names = TRUE)
ae_categories <- purrr::map_dfr(ae_files, function(f) {
  df <- read.table(
    f, sep = "\t", header = FALSE,
    col.names = c("orthogroupID", "ensemblID", "category", "genesymbol"),
    stringsAsFactors = FALSE, fill = TRUE
  )
  
  parts <- sub("_AEcateg.tsv$", "", basename(f))
  species <- strsplit(parts, "_")[[1]][1]
  tissue  <- strsplit(parts, "_")[[1]][2]
  
  df %>%
    mutate(
      species = species,
      tissue = tissue,
      category = gsub("–", "-", category)
    )
}) %>%
  filter(!is.na(category), category != "")

print("AE categories:")
print(table(ae_categories$category))

# 2. Load ALL rewired edges (from updated *_rewiring_all_edges.tsv)
all_rewired_files <- list.files(GRN_DIR, pattern = "_rewiring_all_edges.tsv$", full.names = TRUE)
all_rewired_edges <- map_dfr(all_rewired_files, ~ {
  df <- read.table(
    .x, 
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE   # preserve the ugly SpeciesTissue etc.
  )
  
  # Extract Species and Tissue from the filename
  parts <- strsplit(sub("_rewiring_all_edges.tsv$", "", basename(.x)), "_")[[1]]
  species <- parts[1]
  tissue  <- parts[2]
  
  df %>%
    mutate(
      Tissue = tissue,
      Species = species,
      Species_Tissue = paste0(tissue, "_", species)
    )
})

# 3. Annotate ALL edges with 4 categories
hahe_ids <- unique(filter(ae_categories, category == "HA-HE")$orthogroupID)
mame_ids <- unique(filter(ae_categories, category == "MA-ME")$orthogroupID)

all_rewired_edges <- all_rewired_edges %>%
  mutate(
    Reg_cat = case_when(
      Regulator_Orthogroup %in% hahe_ids ~ "HA-HE",
      Regulator_Orthogroup %in% mame_ids ~ "MA-ME", TRUE ~ "Other"
    ),
    Tgt_cat = case_when(
      Target_Orthogroup %in% hahe_ids ~ "HA-HE",
      Target_Orthogroup %in% mame_ids ~ "MA-ME", TRUE ~ "Other"
    ),
    Edge_cat = case_when(
      Reg_cat == "HA-HE" & Tgt_cat == "HA-HE" ~ "HA-HE → HA-HE",
      Reg_cat == "HA-HE" & Tgt_cat != "HA-HE" ~ "HA-HE Regulator",
      Reg_cat != "HA-HE" & Tgt_cat == "HA-HE" ~ "HA-HE Target",
      Reg_cat == "MA-ME" & Tgt_cat == "MA-ME" ~ "MA-ME only",
      TRUE ~ "Other"
    ),
    Edge_cat = factor(Edge_cat, levels = c("HA-HE → HA-HE", "HA-HE Regulator", "HA-HE Target", "MA-ME only", "Other"))
  )

print("ALL rewired edges loaded:")
print(all_rewired_edges %>% 
        count(Species, Tissue, wt = NULL, name = "n_edges"))

# 4. Create top 12k "foreground" (highest rewiring_score within ALL)
top12k_rewired <- all_rewired_edges %>%
  group_by(Species_Tissue) %>%
  slice_max(rewiring_score, n = TOP_REWIRED, with_ties = FALSE) %>%  # Top 12k by rewiring
  ungroup()

print("Top 12k rewired counts:")
print(top12k_rewired %>% count(Species, Tissue, Species_Tissue))

print("Category counts (ALL vs top 12k):")
print(bind_rows(
  all_rewired_edges %>% count(Edge_cat, name = "all_n"),
  top12k_rewired %>% count(Edge_cat, name = "top12k_n")
))

# Edge_cat  all_n top12k_n
# 1    HA-HE → HA-HE 866900       NA
# 2  HA-HE Regulator 603445       NA
# 3     HA-HE Target 801321       NA
# 4       MA-ME only 124569       NA
# 5            Other 596587       NA
# 6    HA-HE → HA-HE     NA   102875
# 7  HA-HE Regulator     NA    47828
# 8     HA-HE Target     NA    49932
# 9       MA-ME only     NA     5784
# 10           Other     NA    21581

# 5. Enrichment test: top 12k rewired (FG) vs rest of ALL rewired (BG)
all_edge_levels <- c("HA-HE → HA-HE", "HA-HE Regulator", "HA-HE Target", "MA-ME only")

enrichment_results <- map_dfr(unique(all_rewired_edges$Species_Tissue), ~ {
  fg_df <- filter(top12k_rewired, Species_Tissue == .x)
  bg_df <- filter(all_rewired_edges, Species_Tissue == .x)
  total_fg <- nrow(fg_df); total_bg_rest <- nrow(bg_df) - total_fg
  
  map2_dfr(.x, all_edge_levels, ~ {
    n_fg_cat <- sum(fg_df$Edge_cat == .y, na.rm = TRUE)
    n_bg_cat <- sum(bg_df$Edge_cat == .y, na.rm = TRUE)
    n_bg_rest_cat <- max(0, n_bg_cat - n_fg_cat)
    n_fg_other <- total_fg - n_fg_cat
    n_bg_rest_other <- total_bg_rest - n_bg_rest_cat
    
    mat <- matrix(c(n_fg_cat, n_fg_other, n_bg_rest_cat, n_bg_rest_other), nrow = 2)
    
    if (sum(mat) == 0 || min(mat) < 1) {
      ft_p <- 1; ft_or <- 1
    } else {
      ft <- fisher.test(mat, alternative = "greater")
      ft_p <- ft$p.value; ft_or <- unname(ft$estimate)
    }
    
    tibble(
      Species_Tissue = .x,  # Species_Tissue from outer loop
      Edge_cat = .y,        # Edge level from all_edge_levels
      total_top12k_rewired = total_fg, total_rest_all_rewired = total_bg_rest,
      n_top12k_rewired = n_fg_cat, n_rest = n_bg_rest_cat,
      prop_top12k_rewired = n_fg_cat / total_fg,
      prop_all_rewired = n_bg_cat / nrow(bg_df),
      enrichment = (n_fg_cat / total_fg) / (n_bg_cat / nrow(bg_df)),
      fisher_p = ft_p, odds_ratio = ft_or
    )
  })
}) %>%
  group_by(Species_Tissue) %>%
  mutate(fisher_p_adj = p.adjust(fisher_p, method = "BH")) %>%
  ungroup() %>%
  arrange(Species_Tissue, Edge_cat, desc(enrichment))

write_tsv(enrichment_results, file.path(AE_DIR, "HAHE_enrichment_top12k_rewired_vs_all_rewired_4cats.tsv"))
print("Enrichment (top 12k rewired vs ALL rewired):")
print(enrichment_results %>% filter(Edge_cat != "Other") %>% select(Species_Tissue, Edge_cat, enrichment, fisher_p, fisher_p_adj))

# 5a. plot the enrichment
library(tidyverse)
library(ggdist)
library(ggrepel)

# 5b. Prepare a long “ranked” version for plotting
p_sig <- 0.05

plot_data <- enrichment_results %>%
  filter(!is.na(enrichment), enrichment > 0) %>%
  arrange(Species_Tissue, desc(enrichment)) %>%
  group_by(Species_Tissue) %>%
  mutate(
    rank_cat = row_number(),  # 1 = most enriched category per sp‑tis
    is_sig = fisher_p_adj < p_sig,
    label = ifelse(
      is_sig,
      paste0(Edge_cat, " (adj‑p = ", scales::number(fisher_p_adj, digits = 2), ")"),
      Edge_cat
    )
  ) %>%
  ungroup()

# 5c. Plot enrichment rank per species‑tissue

library(tidyverse)
library(ggrepel)
library(patchwork)

p_sig <- 0.05

plot_data <- enrichment_results %>%
  mutate(
    enrichment = as.numeric(enrichment),
    fisher_p = as.numeric(fisher_p),
    fisher_p_adj = as.numeric(fisher_p_adj),
    neglog10_padj = -log10(fisher_p_adj),
    is_sig = !is.na(fisher_p_adj) & fisher_p_adj < p_sig
  ) %>%
  filter(!is.na(enrichment), !is.na(fisher_p_adj), enrichment > 0) %>%
  arrange(Species_Tissue, desc(enrichment)) %>%
  group_by(Species_Tissue) %>%
  mutate(
    rank_cat = row_number(),
    label = paste0(Edge_cat, "\n", "adj-p=", signif(fisher_p_adj, 3))
  ) %>%
  ungroup()

# Check that significance is being assigned correctly
print(plot_data %>% select(Species_Tissue, Edge_cat, enrichment, fisher_p_adj, neglog10_padj, is_sig))

# Panel 1: enrichment
p_enrich <- plot_data %>%
  ggplot(aes(x = rank_cat, y = enrichment, color = is_sig, shape = Species_Tissue)) +
  geom_point(size = 2, alpha = 0.9) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  scale_color_manual(
    values = c("TRUE" = "#D55E00", "FALSE" = "#999999"),
    breaks = c(TRUE, FALSE),
    labels = c("Significant (adj-p < 0.05)", "Not significant")
  ) +
  facet_wrap(~ Species_Tissue, scales = "free_x") +
  labs(
    title = "Ranked enrichment of edge categories in top 12k rewired edges",
    subtitle = "Colored by BH-adjusted Fisher test significance",
    x = "Rank within species-tissue",
    y = "Enrichment",
    color = "Significance",
    shape = "Accessibility-Expression Category"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# Panel 2: -log10 adjusted p-value
p_padj <- plot_data %>%
  ggplot(aes(x = rank_cat, y = neglog10_padj, color = is_sig, shape = Species_Tissue)) +
  geom_point(size = 2, alpha = 0.9) +
  geom_hline(yintercept = -log10(p_sig), linetype = "dashed", color = "grey50") +
  scale_color_manual(
    values = c("TRUE" = "#D55E00", "FALSE" = "#999999"),
    breaks = c(TRUE, FALSE),
    labels = c("Significant (adj-p < 0.05)", "Not significant")
  ) +
  facet_wrap(~ Species_Tissue, scales = "free_x") +
  labs(
    x = "Rank within species-tissue",
    y = expression(-log[10](adj~p)),
    color = "Significance",
    shape = "Accessibility-Expression Category"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

p_final <- p_enrich / p_padj + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  file.path(OUT_DIR, "HAHE_ranked_enrichment_by_category_with_padj.png"),
  p_final, width = 16, height = 10, dpi = 300, bg = "white"
)


# 6. Score distributions (half-eye plots on ALL rewired edges)
cbPalette <- c("HA-HE → HA-HE" = "#D55E00", "HA-HE Regulator" = "#E69F00", 
               "HA-HE Target" = "#009E73", "MA-ME only" = "#56B4E9", "Other" = "#999999")

# Per tissue, facet by species
# This gives one figure per measure, with tissues on the x-axis and species as facets:

species_names <- c(
  "Ab" = "A. burtoni",
  "Mz" = "M. zebra",
  "Nb" = "N. brichardi",
  "On" = "O. niloticus",
  "Pn" = "P. nyererei"
)

tissue_names <- c(
  "B" = "Brain",
  "E" = "Retina",
  "L" = "Liver",
  "T" = "Testis"
)

cbPalette <- c(
  "HA-HE → HA-HE" = "#D55E00",
  "HA-HE Regulator" = "#E69F00",
  "HA-HE Target" = "#009E73",
  "MA-ME only" = "#56B4E9",
  "Other" = "#999999"
)

for (tis in sort(unique(all_rewired_edges$Tissue))) {
  
  df_tis <- all_rewired_edges %>%
    filter(Tissue == tis) %>%
    mutate(Species_full = factor(species_names[Species], levels = species_names))
  
  p1_tis <- df_tis %>%
    ggplot(aes(x = Edge_cat, y = CompositeWeight, fill = Edge_cat)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.75, color = NA) +
    geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.35, color = "black") +
    scale_fill_manual(values = cbPalette) +
    facet_wrap(~ Species_full) +
    scale_y_continuous(limits = c(0, 0.06)) +
    labs(
      x = "",
      y = "Composite Weight",
      title = paste0(tissue_names[[tis]], ": ALL Rewired Edges")
    ) +
    theme_bw(base_size = 16) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      legend.position = "none",
      panel.grid.major.x = element_blank()
    )
  
  ggsave(
    filename = file.path(AE_DIR, paste0("HAHE_composite_violin_box_", tis, ".png")),
    plot = p1_tis,
    width = 14,
    height = 8,
    dpi = 300,
    bg = "white"
  )
  
  df_tis2 <- all_rewired_edges %>%
    filter(Tissue == tis, rewiring_score > 0) %>%
    mutate(Species_full = factor(species_names[Species], levels = species_names))
  
  p2_tis <- df_tis2 %>%
    ggplot(aes(x = Edge_cat, y = rewiring_score, fill = Edge_cat)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.75, color = NA) +
    geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.35, color = "black") +
    scale_fill_manual(values = cbPalette) +
    facet_wrap(~ Species_full) +
    scale_y_continuous(limits = c(0, 15)) +
    labs(
      x = "",
      y = "Rewiring Score",
      title = paste0(tissue_names[[tis]], ": ALL Rewired Edges")
    ) +
    theme_bw(base_size = 16) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      legend.position = "none",
      panel.grid.major.x = element_blank()
    )
  
  ggsave(
    filename = file.path(AE_DIR, paste0("HAHE_rewiring_violin_box_", tis, ".png")),
    plot = p2_tis,
    width = 14,
    height = 8,
    dpi = 300,
    bg = "white"
  )
}

# Per tissue, pooled species
# This combines all species within each tissue into a single facet per tissue:

library(tidyverse)
library(patchwork)

species_names <- c(
  "Ab" = "A. burtoni",
  "Mz" = "M. zebra",
  "Nb" = "N. brichardi",
  "On" = "O. niloticus",
  "Pn" = "P. nyererei"
)

tissue_names <- c(
  "B" = "Brain",
  "E" = "Retina",
  "L" = "Liver",
  "T" = "Testis"
)

cbPalette <- c(
  "HA-HE → HA-HE" = "#D55E00",
  "HA-HE Regulator" = "#E69F00",
  "HA-HE Target" = "#009E73",
  "MA-ME only" = "#56B4E9",
  "Other" = "#999999"
)

for (tis in sort(unique(all_rewired_edges$Tissue))) {
  
  df_tis <- all_rewired_edges %>%
    filter(Tissue == tis) %>%
    mutate(
      Species_full = factor(species_names[Species], levels = species_names)
    )
  
  p1 <- df_tis %>%
    ggplot(aes(x = Edge_cat, y = CompositeWeight, fill = Edge_cat)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.75, color = NA) +
    geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.35, color = "black") +
    scale_fill_manual(values = cbPalette) +
    facet_wrap(~ Species_full) +
    scale_y_continuous(limits = c(0, 0.06)) +
    labs(
      x = "",
      y = "Composite Weight",
      title = paste0(tissue_names[[tis]], ": ALL Rewired Edges")
    ) +
    theme_bw(base_size = 16) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      legend.position = "none",
      panel.grid.major.x = element_blank()
    )
  
  p2 <- df_tis %>%
    filter(rewiring_score > 0) %>%
    ggplot(aes(x = Edge_cat, y = rewiring_score, fill = Edge_cat)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.75, color = NA) +
    geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.35, color = "black") +
    scale_fill_manual(values = cbPalette) +
    facet_wrap(~ Species_full) +
    scale_y_continuous(limits = c(0, 15)) +
    labs(
      x = "",
      y = "Rewiring Score",
      title = paste0(tissue_names[[tis]], ": ALL Rewired Edges")
    ) +
    theme_bw(base_size = 16) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      legend.position = "none",
      panel.grid.major.x = element_blank()
    )
  
  p_final <- p1 / p2 + plot_layout(guides = "collect")
  
  ggsave(
    filename = file.path(AE_DIR, paste0("HAHE_", tissue_names[[tis]], "_violin_box_combined.png")),
    plot = p_final,
    width = 14,
    height = 12,
    dpi = 300,
    bg = "white"
  )
}

# 7. Top 50 HA-HE → HA-HE candidates (from top 12k rewired)
candidates <- top12k_rewired %>%
  filter(Edge_cat == "HA-HE → HA-HE") %>%
  mutate(candidate_score = CompositeWeight * pmax(rewiring_score, 0)) %>%
  arrange(desc(candidate_score)) %>%
  slice_head(n = 50) %>%
  select(Species, Tissue, Species_Tissue, Regulator_Orthogroup, Target_Orthogroup,
         Regulator_UnifiedGS, Target_UnifiedGS, CompositeWeight, rewiring_score, 
         candidate_score, everything())

write_tsv(candidates, file.path(AE_DIR, "top50_HAHE_candidates_top12k_rewired.tsv"))

# 8. Summary table
summary_stats <- all_rewired_edges %>%
  group_by(Species_Tissue, Edge_cat) %>%
  summarise(
    n_all = n(),
    mean_comp_all = mean(CompositeWeight, na.rm = TRUE),
    mean_rewire_all = mean(rewiring_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    top12k_rewired %>% 
      count(Species_Tissue, Edge_cat, name = "n_top12k") %>%
      mutate(n_top12k = replace_na(n_top12k, 0)),
    by = c("Species_Tissue", "Edge_cat")
  ) %>%
  mutate(prop_top12k = n_top12k / sum(n_top12k, na.rm = TRUE))  # Within sp-tis

write_tsv(summary_stats, file.path(AE_DIR, "HAHE_summary_all_vs_top12k_rewired.tsv"))

print("✅ COMPLETE: Enrichment top 12k rewired vs ALL rewired + plots + candidates!")

# 9. Figure to show ranking of species-tissue categories by the composite ranking score, with clear visual evidence that the top 30 are all HA-HE associated edges.

library(tidyverse)
library(readxl)
library(scales)
library(janitor)

df <- read_excel(
  "/Users/tarangmehta/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/TGAC/Technology_Development/ATAC-Seq/ATAC_manuscript/Files/genie3/genie3_composite/comparative_GRN/AEanalysis/HAHE_summary_all_vs_top12k_rewired.xlsx",
  sheet = "HAHE_summary_all_vs_top12k_rewi"
) |>
  clean_names() |>
  rename(
    ranking_score = `ranking_score_mean_comp_all_mean_rewire_all_prop_top12k`
  ) |>
  mutate(
    tissue_expanded = case_when(
      str_detect(tissue_species, "^B") ~ "Brain",
      str_detect(tissue_species, "^E") ~ "Retina",
      str_detect(tissue_species, "^L") ~ "Liver",
      str_detect(tissue_species, "^T") ~ "Testis"
    ),
    species_expanded = case_when(
      str_detect(tissue_species, "Ab$") ~ "A. burtoni",
      str_detect(tissue_species, "Mz$") ~ "M. zebra",
      str_detect(tissue_species, "Pn$") ~ "P. nyererei",
      str_detect(tissue_species, "Nb$") ~ "N. brichardi",
      str_detect(tissue_species, "On$") ~ "O. niloticus"
    ),
    full_label = paste0(species_expanded, " ", tissue_expanded),
    rank = dense_rank(desc(ranking_score)),
    category = case_when(
      edge_cat == "HA-HE -> HA-HE" ~ "HA-HE→HA-HE",
      edge_cat == "HA-HE Regulator" ~ "HA-HE Regulator",
      edge_cat == "HA-HE Target" ~ "HA-HE Target",
      TRUE ~ "Other"
    ),
    category = factor(
      category,
      levels = c("HA-HE→HA-HE", "HA-HE Regulator", "HA-HE Target", "Other")
    ),
    top30 = rank <= 30
  )

print(paste("Top 30 all HA-HE:", all(df$top30[1:30] & str_detect(df$edge_cat[1:30], "HA-HE"))))
print(paste("Total categories:", nrow(df)))

p <- ggplot(df, aes(x = reorder(full_label, ranking_score), y = ranking_score, fill = category)) +
  geom_col(width = 0.8, color = "white", linewidth = 0.25) +
  # geom_hline(yintercept = df$ranking_score[30], linetype = "dashed", color = "grey40", linewidth = 0.6) +
  # geom_text(
  #   data = df[df$rank <= 35, ],
  #   aes(label = rank, y = ranking_score + max(ranking_score)*0.02),
  #   size = 4.5, fontface = "bold", vjust = 0
  # ) +
  # annotate("text", x = 1, y = df$ranking_score[30]*1.1,
  #          label = "Top 30\n(all HA-HE)", size = 5, hjust = 0, 
  #          fontface = "bold", color = "#009E73") +
  scale_fill_manual(
    values = c(
      "HA-HE→HA-HE" = "#D55E00",
      "HA-HE Regulator" = "#E69F00",
      "HA-HE Target" = "#009E73",
      "Other" = "#ffe6e6"
    ),
    name = "Edge category"
  ) +
  scale_y_continuous(
    labels = label_scientific(),
    expand = expansion(mult = c(0, 0.15))
  ) +
  labs(
    x = "Species-tissue categories\n(ranked by composite score)",
    y = "Composite Score\n(Mean Composite Weight × Mean Rewiring Score × Proportion of Top 12k rewired edges)"
  ) +
  coord_flip(clip = "off") +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 14),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "bottom",
    plot.margin = margin(10, 10, 10, 10)
  )

# Export for manuscript
ggsave(
  filename = "/Users/tarangmehta/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/TGAC/Technology_Development/ATAC-Seq/ATAC_manuscript/Figures/fig6_or_supp/6c_HAHE_ranking.png",
  plot = p,
  width = 30, height = 12, units = "cm", dpi = 300, bg = "white"
)


###########################
## Fig. 6 Multiplot
###########################

library(cowplot)
library(ggplot2)
library(magick)

p1 <- ggdraw() + draw_image("/Users/tarangmehta/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/TGAC/Technology_Development/ATAC-Seq/ATAC_manuscript/Figures/fig6_or_supp/6a_genie3_rewiringcategories.png")
p2 <- ggdraw() + draw_image("/Users/tarangmehta/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/TGAC/Technology_Development/ATAC-Seq/ATAC_manuscript/Figures/fig6_or_supp/6b_motifGENIE3_rewiringcategories.png")
p3 <- ggdraw() + draw_image("/Users/tarangmehta/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/TGAC/Technology_Development/ATAC-Seq/ATAC_manuscript/Figures/fig6_or_supp/6c_HAHE_ranking.png")

fig6_combined_plot <- plot_grid(
  p1, p2, p3,
  ncol = 1,
  labels = c("A", "B", "C"),
  label_size = 20,
  label_fontface = "bold",
  align = "v"
)

ggsave(
  filename = "/Users/tarangmehta/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/TGAC/Technology_Development/ATAC-Seq/ATAC_manuscript/Figures/fig6_or_supp/fig6_GRNs.png",
  plot = fig6_combined_plot,
  width = 18, height = 26, units = "cm", dpi = 300, bg = "white"
)
