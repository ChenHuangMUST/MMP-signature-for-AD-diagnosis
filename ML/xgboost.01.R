load("D:/RProject/AD/AD_train.RData")
load("D:/RProject/AD/NMF/afterbatch.marker.select.3.11")
load("D:/RProject/AD/AD-reset/NMF/imp_MPP.RData")
library(xgboost)
library(Matrix)
library(pROC)

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

out<-matrix(nrow = 1000,ncol = 4000)
for (i in 1:1000){
  set.seed(i)
ind <- sample(2, nrow(inner), replace = TRUE, prob = c(0.7, 0.3))
inner[,1] <- factor(inner[,1])
train <- inner[ind==1, ]
test <- inner[ind==2, ]
train1=train

#将自变量转化为矩阵
traindata1 <- data.matrix(train[,c(2:ncol(train))])

#利用Matrix函数，将sparse参数设置为TRUE，转化为稀疏矩阵
traindata2 <- Matrix(traindata1,sparse = T)
#将因变量转化为numeric
traindata3 <- as.numeric(as.character(train[,1]))
#将自变量和因变量拼接为list
traindata4 <- list(data=traindata2,label=traindata3)
#构造模型需要的xgd.DMatrix对象，处理对象为稀疏矩阵
dtrain <- xgb.DMatrix(data=traindata4$data,label=traindata4$label)

#定义模型参数
param <- list(max_depth=2,eta=1,silent=1,objective='binary:logistic')
#构造xgboost模型
bst=xgb.train(params = param,data=dtrain,nrounds = 4,nthread=2)
#显示计算过程，查看树特征
model <- xgb.dump(bst,with_stats = T)

#获取特征的真实名称
names <- dimnames(data.matrix(train[,c(2:ncol(train))]))[2]
#计算变量重要性
importance_matrix <- xgb.importance(feature_names=colnames(bst$feature_names),model = bst)


traindata1 <- data.matrix(train[,c(2:ncol(test))])

#利用Matrix函数，将sparse参数设置为TRUE，转化为稀疏矩阵
traindata2 <- Matrix(traindata1,sparse = T)
#将因变量转化为numeric
traindata3 <- as.numeric(as.character(train[,1]))
#将自变量和因变量拼接为list
traindata4 <- list(data=traindata2,label=traindata3)
#构造模型需要的xgd.DMatrix对象，处理对象为稀疏矩阵
dtrain <- xgb.DMatrix(data=traindata4$data,label=traindata4$label)
#模型训练
xgb <- xgboost(data = dtrain,max_depth=6, eta=0.5,  
               objective='binary:logistic', nround=25)
importance <-data.frame( xgb.importance(model = xgb))  

testdata1 <- data.matrix(test[,c(2:ncol(test))])

testdata2 <- Matrix(testdata1,sparse = T)
#将因变量转化为numeric
testdata3 <- as.numeric(as.character(test[,1]))
#将自变量和因变量拼接为list
testdata4 <- list(data=testdata2,label=testdata3)
#构造模型需要的xgd.DMatrix对象，处理对象为稀疏矩阵
dtest <- xgb.DMatrix(data=testdata4$data,label=testdata4$label)
#混淆矩阵
pre_tr = round(predict(xgb,newdata = dtrain))
pre_xgb = round(predict(xgb,newdata = dtest))
table(testdata3,pre_xgb,dnn=c("true","pre"))
xgboost_tr
#ROC曲线
xgboost_tr <- roc(traindata3,as.numeric(pre_tr))
xgboost_roc <- roc(testdata3,as.numeric(pre_xgb))

out[i,] <-c(xgboost_tr$auc,xgboost_roc$auc,importance[,1],rep(NA,time=3998-nrow(importance)))
}
save(out,file="xgboost.RData")
save.image()
