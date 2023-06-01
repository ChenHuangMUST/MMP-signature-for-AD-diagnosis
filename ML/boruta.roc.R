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
