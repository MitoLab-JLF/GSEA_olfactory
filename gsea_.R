# Libraries loading
library(fgsea)
library(data.table)

# GMT file loading
odor_sets <- gmtPathways("olf_upsit_categories.gmt")

#  Deleting of empty columns
odor_sets <- lapply(odor_sets, function(x) x[x != ""])

# RNK file loading and preparation
rnk_data  <- read.table("olf_odors_RS.rnk", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

odor_ranks <- rnk_data$V2
names(odor_ranks) <- rnk_data$V1

# Sorting
odor_ranks <- sort(odor_ranks, decreasing = FALSE)

# --- REPRODUCIBILITY ---
set.seed(42)

# GSEA analysis
fgseaRes <- fgsea(pathways = odor_sets, 
                  stats    = odor_ranks,
                  minSize  = 3, 
                  maxSize  = 50,)

# Results
View(fgseaRes)