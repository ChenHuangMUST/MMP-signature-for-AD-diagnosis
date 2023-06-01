#LASSO####
load("D:/RProject/AD/AD_train.RData")
load("D:/RProject/AD/NMF/afterbatch.marker.select.3.11")
load("D:/RProject/AD/AD-reset/NMF/imp_MPP.RData")
index=imp_big
dat=dat_rmbatch[index[,1],]
#dat_nondense=dat
#dat_nondense[dat_nondense<0] <- -1
#dat_nondense[dat_nondense>0] <- 1
index=intersect(pd[,1],colnames(dat_rmbatch))
pd=pd[index,]
pathway=dat[,index]
pathway=t(pathway)

dat1=cbind(pd,pathway)
#dat_nondense=cbind(pd,t(dat_nondense))
dat1=dat1[,-1]
dat1=apply(dat1,2,as.numeric)
row.names(dat1)=row.names(pathway)
inner=dat1
inner=data.frame(inner)

colnames(inner)[1]="status"
set.seed(2023414)
ind <- sample(2, nrow(inner), replace = TRUE, prob = c(0.7, 0.3))

inner[,1] <- factor(inner[,1])
train <- inner[ind==1, ]
test <- inner[ind==2, ]
#dat_nondense=dat
#dat_nondense[dat_nondense<0] <- -1
#dat_nondense[dat_nondense>0] <- 1

y=as.matrix(inner[,1])
x=as.matrix(inner[,2:ncol(inner)])
x=apply(x, 2, as.numeric)
row.names(x)=row.names(inner)
row.names(y)=row.names(inner)
library("glmnet")
library("survival")
350
out<-matrix(nrow = 5000,ncol = 250)
for (i in 1:5000) {
  set.seed(i)
  index <- sample(nrow(x),size = round(nrow(x)*0.7))
  df=cbind(as.factor(y),x)
  colnames(df)[1]="status"
  df.training=as.matrix(df[index,])
  df.test=as.matrix(df[-index,])
  X_train <- as.matrix(x[index,])
  y_train=as.matrix(y[index,])
  X_test <- as.matrix(x[-index,])
  y_test=as.matrix(y[-index,])
  
  
  fit_train=glmnet(X_train, y_train, family = "binomial", maxit = 1000)
  plot(fit_train, xvar = "lambda", label = TRUE)
  
  cvfit_train= cv.glmnet(X_train, y_train, family="binomial", maxit = 1000)
  plot(cvfit_train)
  #其中两条虚线分别指示了两个特殊的λ值
  dev.off()
  ###4. 输出预测模型的相关系数与riskScore
  ###4.1 输出相关系数
  coef_train=coef(fit_train, s = cvfit_train$lambda.1se)
  index=which(coef_train != 0)
  index=index[-which(index==1)]
  if (length(index)>1) {
    index=index
  }else{index=c(2,3)}
  actCoef_train=coef_train[index]
  lassoGene_train=row.names(coef_train)[index]
  if(length(lassoGene_train)!=0) {
    lassoGene_train=lassoGene_train
  } else {
    lassoGene_train=c("hsa00531-hsa00640","hsa00230-hsa00240")
  }
  geneCoef_train=cbind(Gene=lassoGene_train,Coef=actCoef_train)
  geneCoef_train   #查看模型的相关系数
  index
  ###4.2 计算riskScore
  FinalGeneExp_train =df.training[,lassoGene_train]
  myFun_train = function(x){crossprod(as.numeric(x),actCoef_train)}
  riskScore_train = apply(FinalGeneExp_train,1,myFun_train)
  outCol_train = c("status", lassoGene_train)
  risk_train = as.vector(ifelse(riskScore_train > median(riskScore_train), "high", "low"))
  c_train=apply(df.training[,outCol_train], 2,as.numeric)
  row.names(c_train)=row.names(df.training)
  ###5. 绘制散点分布图
  #install.packages("ggpubr")
  library(ggpubr)  
  
  c_train=as.data.frame(c_train)
  c_train$risk=risk_train
  c_train$riskScore=riskScore_train
  dat=c_train
  
  
  ###6. 判断预测结果的准确性
  #install.packages("ROCR")
  library(ROCR)   #使用ROCR包绘制预测模型的ROC曲线
  library(glmnet)
  library(caret)
  
  pred <- prediction(dat$riskScore, dat$status)
  perf <- performance(pred,"tpr","fpr")
  AUC <- performance(pred,"auc")   #计算AUC
  plot(perf,colorize=FALSE, col="red", print.auc =TRUE) #绘制ROC曲线
  lines(c(0,1),c(0,1),col = "gray", lty = 4 )
  AUC@y.values
  a1=as.numeric(AUC@y.values)
  FinalGeneExp_test =df.test[,lassoGene_train]
  fit_test=glmnet(X_test, y_test, family = "binomial", maxit = 1000)
  cvfit_test= cv.glmnet(X_test, y_test, family="binomial", maxit = 1000)
  
  coef_test=coef(fit_test, s = cvfit_test$lambda.min)
  index=which(coef_test != 0)
  actCoef_test=coef_test[index]
  
  
  myFun_test = function(x){crossprod(as.numeric(x),actCoef_train)}
  riskScore_test = apply(FinalGeneExp_test,1,myFun_test)
  outCol_test = c("status", lassoGene_train)
  risk_test = as.vector(ifelse(riskScore_test > median(riskScore_test), "high", "low"))
  c_test=apply(df.test[,outCol_test], 2,as.numeric)
  row.names(c_test)=row.names(df.test)
  
  c_test=as.data.frame(c_test)
  c_test$risk=risk_test
  c_test$riskScore=riskScore_test
  dat=c_test
  pred <- prediction(dat$riskScore, dat$status)
  perf <- performance(pred,"tpr","fpr")
  AUC <- performance(pred,"auc")   #计算AUC
  AUC@y.values
  a2=as.numeric(AUC@y.values)
  out[i,] <- c(a1,a2,lassoGene_train,rep(NA,time=250-(length(lassoGene_train)+2)))
}

write.csv(table(out),file="LASSO.csv")

write.csv(table(out[,3:ncol(out)])/5,file="LASSO_5.csv")
save.image()
load(".RData")
#外部验证
load("D:/RProject/AD/out.test.RData")
{

  outset_gsea=t(outset_gsea)
  otest1=data.frame(apply(outset_gsea,2,as.numeric))
  otest1[,1]=matrix(outset_gsea[,1])
  otest1[,1]=as.factor(otest1[,1])
  rownames(otest1)=rownames(outset_gsea)
  colnames(otest1)=gsub("-",".",colnames(otest1))
  df.test=data.frame(otest1)
  FinalGeneExp_test =df.test[,lassoGene_train]
  X_test=as.matrix(df.test[,2:ncol(df.test)])
  y_test=as.matrix(df.test[,-(2:ncol(df.test))])
  fit_test=glmnet(X_test, y_test, family = "binomial", maxit = 1000)
  cvfit_test= cv.glmnet(X_test, y_test, family="binomial", maxit = 1000)
  
  coef_test=coef(fit_test, s = cvfit_test$lambda.min)
  index=which(coef_test != 0)
  actCoef_test=coef_test[index]
  
  
  myFun_test = function(x){crossprod(as.numeric(x),actCoef_train)}
  riskScore_test = apply(FinalGeneExp_test,1,myFun_test)
  outCol_test = c("status", lassoGene_train)
  risk_test = as.vector(ifelse(riskScore_test > median(riskScore_test), "high", "low"))
  c_test=apply(df.test[,outCol_test], 2,as.numeric)
  row.names(c_test)=row.names(df.test)
  
  c_test=as.data.frame(c_test)
  c_test$risk=risk_test
  c_test$riskScore=riskScore_test
  dat=c_test
  pred <- prediction(dat$riskScore, dat$status)
  perf <- performance(pred,"tpr","fpr")
  AUC <- performance(pred,"auc")   #计算AUC
  AUC@y.values
  a3=as.numeric(AUC@y.values)
}
{
  otest2=data.frame(apply(outset_gsea2,2,as.numeric))
  otest2[,1]=outset_gsea2[,1]
  otest2[,1]=as.factor(otest2[,1])
  colnames(otest2)=gsub("-",".",colnames(otest2))
  df.test=data.frame(otest2)
  colnames(df.test)[1]="status"
  FinalGeneExp_test =df.test[,lassoGene_train]
  X_test=as.matrix(df.test[,2:ncol(df.test)])
  y_test=as.matrix(df.test[,1])
  fit_test=glmnet(X_test, y_test, family = "binomial", maxit = 1000)
  cvfit_test= cv.glmnet(X_test, y_test, family="binomial", maxit = 1000)
  
  coef_test=coef(fit_test, s = cvfit_test$lambda.min)
  index=which(coef_test != 0)
  actCoef_test=coef_test[index]
  
  
  myFun_test = function(x){crossprod(as.numeric(x),actCoef_train)}
  riskScore_test = apply(FinalGeneExp_test,1,myFun_test)
  outCol_test = c("status", lassoGene_train)
  risk_test = as.vector(ifelse(riskScore_test > median(riskScore_test), "high", "low"))
  c_test=apply(df.test[,outCol_test], 2,as.numeric)
  row.names(c_test)=row.names(df.test)
  
  c_test=as.data.frame(c_test)
  c_test$risk=risk_test
  c_test$riskScore=riskScore_test
  dat=c_test
  pred <- prediction(dat$riskScore, dat$status)
  perf <- performance(pred,"tpr","fpr")
  AUC <- performance(pred,"auc")   #计算AUC
  AUC@y.values
  a3=as.numeric(AUC@y.values)
}

save(out,file="lasso.roc.RData")

save.image()
