index=c("hsa00100.hsa00190",
  "hsa00563.hsa00190",
  "hsa00534.hsa00190",
  "hsa00900.hsa00190",
  "hsa00310.hsa00534",
  "hsa00760.hsa00190",
  "hsa00531.hsa00860",
  "hsa00513.hsa00620",
  "hsa01040.hsa00190",
  "hsa00310.hsa00600",
  "hsa00534.hsa00620",
  "hsa00310.hsa00531",
  "hsa00051.hsa00860")
library(pROC)
gl_m <- glm(status ~.,data = test, family = binomial(link = "logit"))
pp=predict.glm(gl_m,newdata = train)
rocpp=roc(train$status,pp)
po=predict.glm(gl_m,newdata = )