library(fgsea)

# 1. Data loading and ranking
odors_DB <- gmtPathways("olf_upsit_categories.gmt")
ranks_df <- read.table("olf_odors_RS.rnk", header = FALSE, sep = "\t", col.names = c("Odor", "Score"))
ranks_df <- ranks_df[order(ranks_df$Score, ranks_df$Odor, decreasing = TRUE), ]

pathway_names <- names(odors_DB)
n_simulations <- 5000

# 2. NES and p-adj matrix
nes_matrix  <- matrix(0, nrow = length(pathway_names), ncol = n_simulations, dimnames = list(pathway_names, NULL))
padj_matrix <- matrix(1, nrow = length(pathway_names), ncol = n_simulations, dimnames = list(pathway_names, NULL))

# --- Permutation cyclus ---
cat("Running 5000 simultions. Collecting NES and p-values...\n")
set.seed(42)

for(i in 1:n_simulations) {
  shuffled_names <- sample(ranks_df$Odor)
  shuffled_ranks <- ranks_df$Score
  names(shuffled_ranks) <- shuffled_names
  
  res <- suppressWarnings(
    fgsea(pathways = odors_DB, stats = shuffled_ranks, minSize = 3, maxSize = 50, gseaParam = 2)
  )
  
  if(!is.null(res) && nrow(res) > 0) {
    pos <- match(res$pathway, pathway_names)
    
    # NES
    valid_nes <- !is.na(pos) & !is.na(res$NES)
    if(any(valid_nes)) {
      nes_matrix[pos[valid_nes], i] <- res$NES[valid_nes]
    }
    
    # Adjusted p-values (padj)
    valid_padj <- !is.na(pos) & !is.na(res$padj)
    if(any(valid_padj)) {
      padj_matrix[pos[valid_padj], i] <- res$padj[valid_padj]
    }
  }
}

nes_matrix[is.na(nes_matrix)] <- 0
padj_matrix[is.na(padj_matrix)] <- 1

# --- Final table with adj values ---
adjusted_simulation_output <- data.frame(
  Category = pathway_names,
  
  # statistics NES
  Sim_Mean_NES = rowMeans(nes_matrix),
  Sim_SD_NES   = apply(nes_matrix, 1, sd),
  
  #  p-adj distribution
  Sim_Min_FDR    = apply(padj_matrix, 1, min),
  Sim_Median_FDR = apply(padj_matrix, 1, median),
  Sim_Max_FDR    = apply(padj_matrix, 1, max),
  
  # Percentage of runs where chance passed the FDR filter (FDR < 0.05 or FDR < 0.25)
  Pct_FDR_under_0.05 = rowSums(padj_matrix <= 0.05) / n_simulations * 100,
  Pct_FDR_under_0.25 = rowSums(padj_matrix <= 0.25) / n_simulations * 100
)

adjusted_simulation_output[,-1] <- round(adjusted_simulation_output[,-1], 4)
adjusted_simulation_output <- adjusted_simulation_output[order(adjusted_simulation_output$Category), ]

cat("\nSimultion finished!\n")
print(adjusted_simulation_output, row.names = FALSE)

write.csv(adjusted_simulation_output, "GSEA_Simulation_Adjusted_Metrics.csv", row.names = FALSE)