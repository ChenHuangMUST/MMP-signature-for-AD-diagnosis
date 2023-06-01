load("D:/RProject/AD/AD_train.RData")
load("D:/RProject/AD/NMF/afterbatch.marker.select.3.11")
load("D:/RProject/AD/AD-reset/NMF/imp_MPP.RData")
dat=dat_rmbatch
#dat_nondense=dat
#dat_nondense[dat_nondense<0] <- -1
#dat_nondense[dat_nondense>0] <- 1
index=intersect(pd[,1],colnames(dat_rmbatch))
pd=pd[index,]
pathway=dat[,index]
pathway=t(pathway)
pathway[pathway<0] <- -1
pathway[pathway>=0] <- 1
dat1=cbind(pd,pathway)
#dat_nondense=cbind(pd,t(dat_nondense))
dat1=dat1[,-1]
dat1=apply(dat1,2,as.numeric)
row.names(dat1)=row.names(pathway)
inner=dat1
inner=data.frame(inner)


library(Boruta)
#install.packages("Boruta")

out<-matrix(nrow = 1000,ncol = 250)
for (i in 1:1000) {
  set.seed(i)
  colnames(inner)[1]="status"
  ind <- sample(2, nrow(inner), replace = TRUE, prob = c(0.7, 0.3))
  inner[,1] <- factor(inner[,1])
  train <- inner[ind==1, ]
  test <- inner[ind==2, ]
  boruta <- Boruta(status~.,data = train,doTrace=2)
  table(boruta$finalDecision)
  boruta$finalDecision
  #提取重要的变量和可能重要的变量
  boruta.finalVarsWithTentative <- data.frame(Item=getSelectedAttributes(boruta, withTentative = T), Type="Boruta_with_tentative")
  out[i,] <-c(boruta.finalVarsWithTentative$Item,rep(NA,time=250-length(boruta.finalVarsWithTentative$Item)))
}

write.csv(table(out),file="Boruta.csv")
save.image()
