library(rpart)
library(pROC)
library(caret)
###
inner[,1] <- factor(inner[,1])
###
out_lst<-list()
for (j in 1:10){
  set.seed(j+100)
  out<-matrix(nrow = 10,ncol = 6)
  folds <- createFolds(pd[,2], k = 10)
  
  for (i in 1:10) {
    
    ind<-folds[[i]]
    
    train <- inner[-ind, ]
    test <- inner[ind, ]
    model<-rpart(status ~.,data=train)
    model <- prune(model, cp = model$cptable[which.min(model$cptable[,"xerror"]),"CP"])
    
    pp <- predict(model , newdata=test,type = "prob")
    pp <- ifelse(pp[,2] > pp[,1],1,0)
    roc_te<-roc((as.numeric(test$status)-1), pp)
    #pp <- round(pp)
    mse <- mean((pp - (as.numeric(test$status)-1))^2)
    cm <- confusionMatrix(data = as.factor(pp), reference = test$status)
    accuracy <- cm$overall["Accuracy"]
    precision <- cm$byClass["Pos Pred Value"]
    recall <- cm$byClass["Sensitivity"]
    f1 <- 2 * precision * recall / (precision + recall)
    out[i,] <-c(roc_te$auc,accuracy,precision,recall,f1,mse)
  }
  
  out_lst[[j]]<-out
}
###

