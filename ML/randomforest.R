load("D:/RProject/AD/AD_train.RData")
load("D:/RProject/AD/NMF/afterbatch.marker.select.3.11")
load("D:/RProject/AD/AD-reset/NMF/imp_MPP.RData")
library(randomForest)
library(pROC)
index=imp_MPP[which(as.numeric(imp_MPP[,4])<0.05),]
dat=dat_rmbatch[index[,1],]
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

colnames(inner)[1]="status"
set.seed(2023414)
out<-matrix(nrow = 1000,ncol = 250)
for (i in 1:1000) {
  set.seed(i)
  ind <- sample(2, nrow(inner), replace = TRUE, prob = c(0.7, 0.3))
  inner[,1] <- factor(inner[,1])
  train <- inner[ind==1, ]
  test <- inner[ind==2, ]
  train1=train
  names(train1) <- make.names(names(train1))
  ran <- randomForest(status~.,data = train1,ntree=500,importance=TRUE,mtry=3,proximity=TRUE)
  ran1 <- ran$importance
  ran1 <- as.data.frame(ran1)
  ran2 <- ran1[order(ran1$MeanDecreaseGini,decreasing = T),]
  ran3 <- subset(ran2,ran2$MeanDecreaseGini>0.1)
  test1=test
  names(test1) <- make.names(names(test1))
  pre_ran <- predict(ran,newdata=test1)
  pre_ran1=predict(ran,newdata=train1)
  obs_ran <- data.frame(prob=pre_ran,obs=test1$status)
  table(test1$status,pre_ran,dnn = c("真实值","预测值"))
  ran_roc <- roc(test1$status,as.numeric(pre_ran))
  ran_roc_tr=roc(train1$status,as.numeric(pre_ran1))

  out[i,] <-c(ran_roc_tr$auc,ran_roc$auc,rownames(ran3),rep(NA,time=248-length(rownames(ran3))))
}
save(out,"randomforest.roc.RData")
write.csv(table(out[,2:ncol(out)]),file="Randomforest.csv")
save.image()