#randomforest

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
