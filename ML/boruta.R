
library(Boruta)
#install.packages("Boruta")
inner[,1] <- factor(inner[,1])

  ind <- sample(2, nrow(inner), replace = TRUE, prob = c(0.7, 0.3))
  inner[,1] <- factor(inner[,1])
  train <- inner[ind==1, ]
  test <- inner[ind==2, ]
  boruta <- Boruta(status~.,data = train,doTrace=2)
  table(boruta$finalDecision)
  boruta$finalDecision
  #Extract important variables and potentially important variables
  boruta.finalVarsWithTentative <- data.frame(Item=getSelectedAttributes(boruta, withTentative = T), Type="Boruta_with_tentative")

#Extract the variables and calculate
#model <- glm(status ~ ., data = train, family = binomial)


