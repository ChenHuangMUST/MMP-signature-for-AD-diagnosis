#BiocManager::install("ReactomePA")
group=read.delim("D:/RProject/AD/NMF/2cluster.consensus.csv",sep=",")
load("D:/RProject/AD/NMF/dat.nmf.ad.RData")
#save.image()
load(".RData")
library(ReactomePA)
library(tidyverse)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(limma)
library(ggplot2)
R.utils::setOption("clusterProfiler.download.method",'auto')
setwd("/Users/fengwen/Desktop/GEOdetails")
save(aaa,choose_matrix,contrast.matrix,DEG,design,exp,fit,fit2,
     g,gene_map,gene.df,genelist_input,group,file="go_kegg.Rdata")
load("go_kegg.Rdata")
#limma差异分析流程
exp=combat_edata
exp=exp[,group$X]
dim(exp)

group=t(cbind(t(groupq),t(groupw)))

group$result_rmb..2.....consensusClass...=gsub("1", "group1", group, perl=TRUE)
group$result_rmb..2.....consensusClass...=gsub("2", "group2", group, perl=TRUE)
design <- model.matrix(~0+factor(group))
colnames(design)=levels(factor(group))
rownames(design)=colnames(exp)
design
##构建差异比较矩阵
contrast.matrix ##这个矩阵声明，我们要把progres.组跟stable进行差异分析比较
##
contrast.matrix<- makeContrasts(group1-group2, levels = design)##原始代码有误，这条是修改
##使用limma进行差异分析
##step1
fit <- lmFit(exp,design)
##step2
fit2 <- contrasts.fit(fit, contrast.matrix) ##这一步很重要
fit2 <- eBayes(fit2)  ## default no trend !!!
##eBayes() with trend=TRUE
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput) 
#write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
head(nrDEG)
##进一步看一下热图

## heatmap 
library(pheatmap)
choose_gene=names(tail(sort(apply(exp,1,mad)),200))
choose_matrix=exp[choose_gene,]
choose_matrix=t(scale(t(choose_matrix)))
pheatmap(choose_matrix)
##粗制火山图
plot(nrDEG$logFC, -log10(nrDEG$P.Value))

##好看一点的

## volcano plot
DEG=nrDEG
logFC_cutoff <- with(DEG,mean(abs( logFC)) + 2*sd(abs( logFC)) )
DEG$change = as.factor(ifelse(DEG$P.Value < 0.05 & abs(DEG$logFC) > logFC_cutoff,
                              ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
this_tile <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                    '\nThe number of up gene is ',nrow(DEG[DEG$change =='UP',]) ,
                    '\nThe number of down gene is ',nrow(DEG[DEG$change =='DOWN',])
)
this_tile
head(DEG)
g = ggplot(data=DEG, aes(x=logFC, y=-log10(P.Value), color=change)) +
  geom_point(alpha=0.4, size=1.75) +
  theme_set(theme_set(theme_bw(base_size=20)))+
  xlab("log2 fold change") + ylab("-log10 p-value") +
  ggtitle( this_tile  ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
  scale_colour_manual(values = c('blue','black','red'))  ## corresponding to the levels(res$change)
print(g)
##取前1000个基因进行分析
gene<-head(rownames(nrDEG),1000)


## get the universal genes and sDEG 转化ID，从SYMBOL到ENTREZID
gene.df <- bitr(gene, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)
head(gene.df)

## KEGG pathway analysis
kk <- enrichKEGG(gene         = gene.df$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
head(kk)[,1:6]

genelist_input=rest
genename <- as.character(genelist_input[,1])
gene_map <- AnnotationDbi::select(org.Hs.eg.db, keys=genename,
                                  keytype="SYMBOL",
                                  columns=c("ENTREZID"))
colnames(genelist_input)=c("Gene","logFC")
colnames(gene_map)=c("Gene","ENTRZID")
aaa<-inner_join(gene_map,genelist_input,by = "Gene")
aaa<-aaa[,-1]
aaa<-na.omit(aaa)
aaa$logFC<-sort(aaa$logFC,decreasing = T)
geneList = aaa[,2]
names(geneList) = as.character(aaa[,1])
kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(kk2)[,1:6]
gseaplot(kk2, geneSetID = "hsa04740 ",by="all")

save(aaa,choose_matrix,contrast.matrix,DEG,design,exp,fit,fit2,g,gene_map,gene.df,genelist_input,group,file="go_kegg.Rdata")

##查看结果，了解基因功能（hsa03030在基因复制通路影响很明显）
R.utils::setOption("clusterProfiler.download.method",'auto')
#GSEA分析——GO
Go_gseresult <- gseGO(geneList, 'org.Hs.eg.db', keyType = "ENTREZID", ont="all", nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
#GSEA分析——KEGG
KEGG_gseresult <- gseKEGG(geneList, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
#GSEA分析——Reactome
Go_Reactomeresult <- gsePathway(geneList, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=0.05)
#保存文件
write.table (Go_gseresult, file ="Go_gseresult.csv", sep =",", row.names =TRUE)
write.table (KEGG_gseresult, file ="KEGG_gseresult.csv", sep =",", row.names =TRUE)
write.table (Go_Reactomeresult, file ="Go_Reactomeresult.csv", sep =",", row.names =TRUE)
Go_gseresult=read.table("Go_gseresult.csv", sep =",")
KEGG_gseresult=read.table("KEGG_gseresult.csv", sep =",")
gseaplot(KEGG_gseresult,1,pvalue_table = TRUE)

#########################################
###        GO Pathway Plot           ### 
##           2021.11.24              ##
######################################

setwd("D:/Your/Working/Directory/")
getwd()

dat = read.table("GO.txt",header = T,sep = "\t")
dat=look
library(ggplot2)#没有自己安装 install.package("ggplot2")
p <- ggplot(dat,aes(y=Gene.ratio,x=Term,fill=PValue)) + 
  geom_bar(stat="identity",position = "dodge") +
  facet_grid(Category~.,scales = "free",space = "free") + 
  coord_flip() + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        strip.text.y = element_text(size = 14),
        legend.position="right",
        legend.title = element_text(size=18),
        legend.text = element_text(size=14),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=18),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14))
p
ggsave(p,filename = "GO.pdf",width = 10,height = 7,dpi=300)
ggsave(p,filename = "GO.jpg",width = 10,height = 7,dpi=300)


