library(xgboost)
library(pROC)

inner[, 1] <- factor(inner[, 1])
train <- inner[ind == 1, ]
test <- inner[ind == 2, ]

# Prepare data for xgboost
dtrain <- xgb.DMatrix(data = train[, -1], label = train$status)
dtest <- xgb.DMatrix(data = test[, -1], label = test$status)

# Define model parameters
param <- list(max_depth = 2, eta = 1, silent = 1, objective = 'binary:logistic')

# Train xgboost model
bst <- xgb.train(params = param, data = dtrain, nrounds = 4, nthread = 2)

# Get feature importances
importance_matrix <- xgb.importance(feature_names = colnames(dtrain), model = bst)
importance <- data.frame(importance_matrix)

# Make predictions on train and test sets
pre_tr <- round(predict(bst, newdata = dtrain))
pre_xgb <- round(predict(bst, newdata = dtest))

# Evaluate model performance using ROC curve
xgboost_tr <- roc(train$status, as.numeric(pre_tr))
xgboost_roc <- roc(test$status, as.numeric(pre_xgb))
