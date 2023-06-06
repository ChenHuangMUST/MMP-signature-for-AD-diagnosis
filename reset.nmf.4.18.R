# Load required data
load(c("imp_MPP.RData", "D:/RProject/AD/AD_train.RData", 
      "D:/RProject/AD/NMF/afterbatch.marker.select.3.11", "MPP.RData"))

library(NMF)

# Subset and preprocess the data
pd <- pd[pd[, 2] == 1, ]
dat_rmbatch <- dat_rmbatch[rownames(imp_MPP[as.numeric(imp_MPP[, 5]) < 0.01, ]), ]
dat <- dat_rmbatch[intersect(rownames(pd), colnames(dat_rmbatch)), ]

# Perform NMF with different ranks
expression_data <- as.matrix(dat, nrow(dat))
nmf_results <- list()
for (i in 2:10) {
  nmf_results[[i]] <- nmf(expression_data, rank = i, nrun = 10, seed = 111111,
                          method = 'brunet', normalize = TRUE)
}

# Define function to create different types of plots
create_plot <- function(nmf_obj, plot_type) {
  plot_file <- paste0(plot_type, "_rank", nmf_obj$rank, ".pdf")
  pdf(file.path("output", plot_file), width = 7, height = 7, onefile = FALSE)
  switch(plot_type,
         "consensusmap" = consensusmap(nmf_obj, annRow = NA, annCol = NA,
                                       main = "Consensus matrix", info = FALSE),
         "coefmap" = coefmap(nmf_obj, annRow = NA, annCol = NA,
                             main = "Metagene contributions in each sample",
                             info = FALSE),
         "basismap" = basismap(nmf_obj, annRow = NA, annCol = NA,
                               main = "Metagenes", info = FALSE),
         "plot" = plot(nmf_obj))
  dev.off()
}

# Create plots for each NMF result
lapply(nmf_results, create_plot, plot_type = "consensusmap")
lapply(nmf_results, create_plot, plot_type = "coefmap")
lapply(nmf_results, create_plot, plot_type = "basismap")
lapply(nmf_results, create_plot, plot_type = "plot")

# Save clustering results
cluster_labels <- predict(nmf_results[[2]])
cluster_labels <- data.frame(cluster = paste0("Cluster", cluster_labels),
                             sample = rownames(cluster_labels))
cluster_labels <- cluster_labels[order(cluster_labels$group), ]
save(cluster_labels, file = "output/clu2-reset.RData")
