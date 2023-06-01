#决策树
load("D:/RProject/AD/AD_train.RData")
load("D:/RProject/AD/NMF/afterbatch.marker.select.3.11")
#dat_nondense=dat
#dat_nondense[dat_nondense<0] <- -1
#dat_nondense[dat_nondense>0] <- 1
index=intersect(pd[,1],colnames(dat_rmbatch))
pd=pd[index,]
pathway=dat_rmbatch[,index]
pathway=t(pathway)
dat1=cbind(pd,pathway)
#dat_nondense=cbind(pd,t(dat_nondense))
dat1=dat1[,-1]
dat1=apply(dat1,2,as.numeric)
row.names(dat1)=row.names(pathway)
inner=dat1
inner=data.frame(inner)

colnames(inner)[1]="status"

out<-matrix(nrow = 1000,ncol = 4000)
for (i in 101:1000) {
  set.seed(i)
  ind <- sample(2, nrow(inner), replace = TRUE, prob = c(0.7, 0.3))
  inner[,1] <- factor(inner[,1])
  train <- inner[ind==1, ]
  test <- inner[ind==2, ]
  model <- train(status ~ ., data=train, method="rpart")
  importance <- varImp(model, scale=FALSE)
  importance=importance$importance
  importance=as.matrix(importance)
  look=importance[which(importance[,1]>0),]
  look=data.frame(look)
  out[i,] <-c(rownames(look),rep(NA,time=4000-length(rownames(look))))
}


write.csv(table(out),file="决策树.csv")
save.image()

###计算roc
library(pROC)
out_roc=matrix(nrow = 1000,ncol = 2)

for (i in 1:1000) {
  index=na.omit(out[i,])
  btrain=train[,index]
  btrain=cbind(train[,1],btrain)
  colnames(btrain)[1]="status"
  btest=test[,index]
  btest=cbind(test[,1],btest)
  colnames(btest)[1]="status"
  gl_m <- glm(status ~.,data = btrain, family = binomial(link = "logit"))
  bb=predict.glm(gl_m,newdata = btrain)
  pp=predict.glm(gl_m,newdata = btest)
  
  roc_te=roc(btest$status,pp)
  roc_tr=roc(btrain$status,bb)
  out_roc[i,]=c(roc_tr$auc,roc_te$auc)
}

out_tl=cbind(out_roc,out)
save.image()