library(lmtest)
library(ggplot2)
library(dplyr)
library(glmnet)

sex_likelihoodratiotest <- lrtest(fit_withsex_relaxed,fit_withoutsex_relaxed )
###no relaxed or glmnet objects can be used with lrtest 

null_test <- lm(total_trainingdata_geneExpression~1, data = total_trainingdata_geneExpression)
lm(pesa_score_responseVariable~1, data = total_trainingdata_geneExpression)



###forward stepwise regression --> used to pick significant features
#the empty model would be a dataset with no coefficients (just pesa score)
##the empty model can be either a numeric or character vector
## data must contain BOTH dependent (pesa score) AND independent (gene ID, sex, age) 



total_trainingdata_geneExpression$pesa_score <- pesa_score_responseVariable
total_trainingdata_geneExpression <- total_trainingdata_geneExpression %>% select(pesa_score, everything())
total_trainingdata_geneExpression$age <- total_trainingdata$Age


empty_forward_model <- lm(pesa_score ~ 1, data = total_trainingdata_geneExpression) 
all(is.na(total_trainingdata_geneExpression$pesa_score))

