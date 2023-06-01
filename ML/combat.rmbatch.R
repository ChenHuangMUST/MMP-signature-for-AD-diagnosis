library(FactoMineR)##没有请先安装
library(factoextra)
library(sva)
library(bladderbatch)
#BiocManager::install("sva")
library(sva)

dat_3AD=apply(dd_3AD,2,as.numeric)
row.names(dat_3AD)=row.names(dd_3AD)
#去批次
combat_edata <- ComBat(dat = dd_3AD, batch = batch_3AD[,2])

#生成pca
pre.pca <- PCA(t(dat_3AD),graph = FALSE)
pre.pca.n <- PCA(t(combat_edata),graph = FALSE)

#绘图
batch_3AD[,2]=c(rep("GSE63060",ncol(dd1)),rep("GSE140829",ncol(dd2)),rep("GSE63061",ncol(dd3)))
fviz_pca_ind(pre.pca,
             geom= "point",
             col.ind = bc_3[,2],
             addEllipses = TRUE,
             legend.title="Group"  )

fviz_pca_ind(pre.pca.n,
             geom= "point",
             col.ind =bc_3[,2],
             addEllipses = TRUE,
             legend.title="Group"  )
