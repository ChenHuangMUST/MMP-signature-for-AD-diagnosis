#load("D:/RProject/AD/heatmap/after.batch.RData")
#load("D:/RProject/AD/AD_train.RData")
#save(combat_edata,pd,file = "af.dat.3.14.rds")
#load(file = "af.dat.3.14.rds")
#save.image("go.kegg.gsea.ad.RData")

load("go.kegg.gsea.ad.RData")
group=read.csv("D:/RProject/AD/NMF/2cluster.consensus.csv")
exp=combat_edata[,group$X]
colnames(group)[2]="group"
row.names(group)=group[,1]
group=group[,-1]
group=gsub("1","group1",group)
group=gsub("2","group2",group)
design <- model.matrix(~0+factor(group))
colnames(design)=levels(factor(group))
rownames(design)=colnames(exp)
 ##这个矩阵声明，我们要把progres.组跟stable进行差异分析比较
##
contrast.matrix<- makeContrasts(group1-group2, levels = design)##原始代码有误，这条是修改
##使用limma进行差异分析
contrast.matrix
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
library(pheatmap)
choose_gene=names(tail(sort(apply(exp,1,mad)),200))
choose_matrix=exp[choose_gene,]
choose_matrix=t(scale(t(choose_matrix)))
pheatmap(choose_matrix)
##粗制火山图
plot(nrDEG$logFC, -log10(nrDEG$P.Value))
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
  geom_text_repel(
    data = DEG[DEG$P.Value<0.05&abs(DEG$logFC)>0.178,],
    aes(label = gene),
    size = 4.5,
    color = "black",
    segment.color = "black", show.legend = FALSE )+
  ggtitle( this_tile  ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
  scale_colour_manual(values = c('blue','black','red'))  ## corresponding to the levels(res$change)
print(g)
##取前1000个基因进行分析
gene=nrDEG[which(nrDEG$logFC>logFC_cutoff&nrDEG$adj.P.Val<0.05),]
library(clusterProfiler)
library(org.Hs.eg.db)
erich.go.BP = enrichGO(gene = row.names(gene),
                       
                       OrgDb = org.Hs.eg.db,
                       
                       keyType = "SYMBOL",
                       
                       pAdjustMethod = "BH",
                       
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)
barplot(erich.go.BP) 
ALL <- enrichGO(gene=row.names(gene),
                
                OrgDb=org.Hs.eg.db,
                
                keyType = "SYMBOL",
                
                ont = 'ALL',
                
                pAdjustMethod = "BH",
                
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)  #readable = TRUE 这行代码可以使enrez转换为gene_symbol输出

library(ggplot2)
barplot(ALL, split="ONTOLOGY")+ facet_grid(ONTOLOGY~.,scale="free") 

kegmt<-read.gmt("C:/Users/pc/Desktop/biotools/c2.cp.kegg.v2022.1.Hs.symbols.gmt") #读gmt文件
KEGG<-GSEA(gene_df,TERM2GENE = kegmt) #GSEA分析
geneList<-gene_df $logFC #第二列可以是folodchange，也可以是logFC
names(geneList)=gene_df $ENTREZID #使用转换好的ID
geneList=sort(geneList,decreasing = T) #从高到低排序
kegmt<-read.gmt("C:/Users/pc/Desktop/biotools/c2.cp.kegg.v2022.1.Hs.symbols.gmt") #读gmt文件
KEGG<-GSEA(gene_df,TERM2GENE = kegmt) #GSEA分析
KEGG<-GSEA(geneList,TERM2GENE = kegmt) #GSEA分析
gene_df <- data.frame(SYMBOL = row.names(nrDEG),
                      logFC=nrDEG$logFC #可以是foldchange
) #记住你的基因表头名字
geneList<-gene_df $logFC #第二列可以是folodchange，也可以是logFC
names(geneList)=gene_df $ENTREZID #使用转换好的ID
geneList=sort(geneList,decreasing = T) #从高到低排序
kegmt<-read.gmt("C:/Users/pc/Desktop/biotools/c2.cp.kegg.v2022.1.Hs.symbols.gmt") #读gmt文件
KEGG<-GSEA(geneList,TERM2GENE = kegmt) #GSEA分析
View(kegmt)
names(geneList)=row.names(nrDEG) #使用转换好的ID
geneList=sort(geneList,decreasing = T) #从高到低排序
kegmt<-read.gmt("C:/Users/pc/Desktop/biotools/c2.cp.kegg.v2022.1.Hs.symbols.gmt") #读gmt文件
KEGG<-GSEA(geneList,TERM2GENE = kegmt) #GSEA分析
library(ggplot2)
dotplot(KEGG) #出点图
dotplot(KEGG,color="pvalue")  #按p值出点图
View(KEGG)
View(KEGG@result)
gogmt=read.gmt("C:/Users/pc/Desktop/biotools/c5.go.bp.v2022.1.Hs.symbols.gmt") #读gmt文件
go<-GSEA(geneList,TERM2GENE = gogmt) #GSEA分析
dotplot(go,color="pvalue")
barplot(go,color="pvalue")
barplot(go)
hmgmt=read.gmt("C:/Users/pc/Desktop/biotools/h.all.v2023.1.Hs.symbols.gmt") #读gmt文件
hm<-GSEA(geneList,TERM2GENE = hmgmt) #GSEA分析
dotplot(hm,color="pvalue")
dotplot(KEGG) #出点图
dotplot(KEGG,color="pvalue")  #按p值出点图
#load("D:/RProject/AD/heatmap/after.batch.RData")
#load("D:/RProject/AD/AD_train.RData")
#save(combat_edata,pd,file = "af.dat.3.14.rds")
#load(file = "af.dat.3.14.rds")
save.image("go.kegg.gsea.ad.RData")
