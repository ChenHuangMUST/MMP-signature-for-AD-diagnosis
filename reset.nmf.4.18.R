load("imp_MPP.RData")
load("D:/RProject/AD/AD_train.RData")
load("D:/RProject/AD/NMF/afterbatch.marker.select.3.11")
load("MPP.RData")
library(NMF)
imp_MPP=imp_big
dat=dat_rmbatch[which(as.numeric(imp_MPP[,5])<0.01),]
index=which(pd[,2]==1)
pd=pd[index,]
index=intersect(rownames(pd),colnames(dat_rmbatch))
dat=dat[,index]
dat_e=exp(1)^dat

x.e <- nmf(dat_e,2:10,nrun=10,seed=111111,method = 'brunet')



plot(x.e)



x.e.nmf.2=nmf(dat_e,2,nrun=10,seed=520,method = 'brunet') # run parameters,rank=3



pdf('consensusmap_x.e.2.1.pdf',width = 7,height = 7,onefile = F)
consensusmap(x.e.nmf.2,
             annRow = NA,
             annCol = NA,
             main = "Consensus matrix",
             info = FALSE)
dev.off() 
pdf('coefmap_rank2.pdf',width = 7,height = 7,onefile = F)
coefmap(x.e.nmf.2,
        annRow = NA,
        annCol = NA,
        main = "Metagene contributions in each sample",
        info = FALSE)
dev.off() 
pdf('basismap_rank2.pdf',width = 7,height = 7,onefile = F)
#作用：每一行显示主导的基底组分，即每一行有最高负载的基底组分。
basismap(x.e.nmf.2,
         annRow = NA,
         annCol = NA,
         main = "Metagenes",
         info = FALSE)
dev.off() 
pdf('plotrank.pdf',width = 7,height = 7,onefile = F)
plot(x.e)
dev.off()

summary(x.e)
group <- predict(x.e.nmf.2)
group <- as.data.frame(group)
group$group <- paste0('Cluster',group$group)
group$sample <- rownames(group)
group<- group[order(group$group),]
table(group$group)
head(group)
save(group,file="clu2-reset.RData")
save.image()
load(".RData")