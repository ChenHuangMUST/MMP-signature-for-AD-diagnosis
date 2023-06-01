load("D:/RProject/AD/NMF/dat.nmf.ad.RData")
library(e1071)
library(preprocessCore)
library(parallel)
library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)
library(ggfortify)
library(GSVA)
kegg_geneset4 <- GSEABase::getGmt("C:/Users/pc/Desktop/biotools/h.all.v2023.1.Hs.symbols.gmt")
hallmarks=read.delim("C:/Users/pc/Desktop/biotools/h.all.v2023.1.Hs.symbols.gmt",header=F)
ssgsea<- gsva(combat_edata,kegg_geneset4,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
cellnum <- ssgsea
cell.prop<- apply(cellnum, 1, function(x){x/sum(x)})
index=intersect(row.names(cell.prop),pd[,1])
cell.prop=cell.prop[index,]
row.names(pd)=pd[,1]
pd=pd[index,]
data4plot <- data.frame()
for (i in 1:ncol(cell.prop)) {
  data4plot <- rbind(
    data4plot,
    cbind(cell.prop[,i],pd[,2],row.names(cell.prop),
          rep(colnames(cell.prop)[i],nrow(cell.prop)
          )
    )
  )
}
data4plot_hallmark=data4plot
identical(row.names(data4plot),row.names(pd))

data4plot
colnames(data4plot)<-c('proportion','group',"sample",'hallmarks')
data4plot$proportion <- as.numeric(data4plot$proportion)


plot_order = data4plot %>% 
  group_by(hallmarks) %>% 
  summarise(m = median(proportion)) %>% 
  arrange(desc(m)) %>% 
  pull(hallmarks)
if(T){
  mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                   axis.title = element_text(size = 12,color ="black"), 
                   axis.text = element_text(size= 12,color = "black"),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   axis.text.x = element_text(angle = 45, hjust = 1 ),
                   panel.grid=element_blank(),
                   legend.position = "top",
                   legend.text = element_text(size= 12),
                   legend.title= element_text(size= 12)
  ) }
box_TME <- ggplot(data4plot, aes(x = hallmarks, y = proportion))+ 
  labs(y="proportion",x= NULL,title = "proportion")+  
  geom_boxplot(aes(fill = group),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+ 
  scale_fill_manual(values = c("#1CB4B8", "#EB7369"))+
  theme_classic() + mytheme + 
  stat_compare_means(aes(group =  group),
                     label = "p.signif",
                     method = "wilcox.test",
                     hide.ns = T)
box_TME
ggplot(data4plot)


library(CIBERSORT)
results <- cibersort(sig_matrix = LM22, mixture_file = combat_edata)

# heatmap
# 按行（样本内部）标准化可以看出在各类样本内部，M2浸润程度（占比）最高
rowscale <- results[,1:ncol(LM22)]#只是相当于备份了一下results
rowscale <- rowscale[,apply(rowscale, 2, function(x){sum(x)>0})]#删除全是0的列
pheatmap(rowscale,
         scale = 'row',#按行标准化，不标准化就会按绝对值显示，很诡异
         cluster_col=T,#是否对列聚类，不聚类，坐标轴就按照原来的顺序显示
         cluster_row=F,#是否对行聚类
         angle_col = "315")#调整X轴坐标的倾斜角度

# 各类样本之间也具有自己占比高的特异性免疫细胞
columnscale <- results[,1:ncol(LM22)]
columnscale <- columnscale[,apply(columnscale, 2, function(x){sum(x)>0})]#删除全是0的列
pheatmap(columnscale,
         scale = 'column',
         cluster_col=F,
         cluster_row=T,
         angle_col = "315")

# 堆积比例图
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87','#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658','#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398','#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963','#968175'
)
cellnum.1 <- results[,1:ncol(LM22)]
cell.prop.1<- apply(cellnum.1, 1, function(x){x/sum(x)})
data4plot.1 <- data.frame()
for (i in 1:ncol(cell.prop)) {
  data4plot.1 <- rbind(
    data4plot.1,
    cbind(cell.prop.1[,i],rownames(cell.prop.1),
          rep(colnames(cell.prop.1)[i],nrow(cell.prop.1)
          )
    )
  )
}


colnames(data4plot.1)<-c('proportion','celltype','sample')
data4plot.1$proportion <- as.numeric(data4plot.1$proportion)
ggplot(data4plot.1,aes(sample,proportion,fill=celltype))+
  geom_bar(stat="identity",position="fill")+
  scale_fill_manual(values=my36colors)+#自定义fill的颜色
  ggtitle("cell portation")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'),axis.title.x=element_text(size=1))+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))+#把x坐标轴横过来
  guides(fill=guide_legend(title=NULL))
iindex=intersect(data4plot.1$sample,row.names(pd))



cell.prop<- apply(results, 1, function(x){x/sum(x)})
cell.prop=t(cell.prop)
index=intersect(row.names(cell.prop),pd[,1])

cell.prop=cell.prop[index,]
row.names(pd)=pd[,1]
pd=pd[index,]
data4plot <- data.frame()
for (i in 1:ncol(cell.prop)) {
  data4plot <- rbind(
    data4plot,
    cbind(cell.prop[,i],pd[,2],row.names(cell.prop),
          rep(colnames(cell.prop)[i],nrow(cell.prop)
          )
    )
  )
}

data4plot <- data.frame()
for (i in 1:ncol(cell.prop)) {
  data4plot <- rbind(
    data4plot,
    cbind(cell.prop[,i],pd[,2],row.names(cell.prop),
          rep(colnames(cell.prop)[i],nrow(cell.prop)
          )
    )
  )
}

identical(row.names(data4plot),row.names(pd))

data4plot
colnames(data4plot)<-c('proportion','group',"sample",'hallmarks')
data4plot$proportion <- as.numeric(data4plot$proportion)


plot_order = data4plot %>% 
  group_by(hallmarks) %>% 
  summarise(m = median(proportion)) %>% 
  arrange(desc(m)) %>% 
  pull(hallmarks)
if(T){
  mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                   axis.title = element_text(size = 12,color ="black"), 
                   axis.text = element_text(size= 12,color = "black"),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   axis.text.x = element_text(angle = 45, hjust = 1 ),
                   panel.grid=element_blank(),
                   legend.position = "top",
                   legend.text = element_text(size= 12),
                   legend.title= element_text(size= 12)
  ) }
box_TME <- ggplot(data4plot, aes(x = hallmarks, y = proportion))+ 
  labs(y="proportion",x= NULL,title = "proportion")+  
  geom_boxplot(aes(fill = group),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+ 
  scale_fill_manual(values = c("#1CB4B8", "#EB7369"))+
  theme_classic() + mytheme + 
  stat_compare_means(aes(group =  group),
                     label = "p.signif",
                     method = "wilcox.test",
                     hide.ns = T)
box_TME
ad