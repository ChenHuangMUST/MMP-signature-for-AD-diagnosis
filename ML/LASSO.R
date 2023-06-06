library("glmnet")
library("survival")

#index was the method of data segmentation

  df <- cbind(as.factor(y), x)
  colnames(df)[1] <- "status"
  df.training <- as.matrix(df[index, ])
  X_train <- as.matrix(x[index, ])
  y_train <- as.matrix(y[index, ])

  fit_train <- glmnet(X_train, y_train, family = "binomial", maxit = 1000)
  cvfit_train <- cv.glmnet(X_train, y_train, family = "binomial", maxit = 1000)

  coef_train <- coef(fit_train, s = cvfit_train$lambda.1se)
  index <- which(coef_train != 0)[-which(which(coef_train != 0) == 1)]
  lassoGene_train <- row.names(coef_train)[index]
  train=df[,lassoGene_train]
  cvfit_train= cv.glmnet(X_train, y_train, family="binomial", maxit = 1000)
