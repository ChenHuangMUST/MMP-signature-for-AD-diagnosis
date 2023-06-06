library(FactoMineR)
library(factoextra)
library(sva)
library(bladderbatch)

# Convert data to numeric format and assign row names
dat_3AD <- apply(dd_3AD, 2, as.numeric)
row.names(dat_3AD) <- row.names(dd_3AD)

# Remove batch effects using ComBat
combat_edata <- ComBat(dat = dd_3AD, batch = batch_3AD[, 2])

# Generate PCA plots for pre-processed and post-processed data
pre.pca <- PCA(t(dat_3AD), graph = FALSE)
pre.pca.n <- PCA(t(combat_edata), graph = FALSE)

# Create group labels for visualization
batch_3AD[, 2] <- c(rep("GSE63060", ncol(dd1)), rep("GSE140829", ncol(dd2)), rep("GSE63061", ncol(dd3)))

# Plot PCA for pre-processed data
fviz_pca_ind(pre.pca,
             geom = "point",
             col.ind = bc_3[, 2],
             addEllipses = TRUE,
             legend.title = "Group")

# Plot PCA for post-processed data
fviz_pca_ind(pre.pca.n,
             geom = "point",
             col.ind = bc_3[, 2],
             addEllipses = TRUE,
             legend.title = "Group")
