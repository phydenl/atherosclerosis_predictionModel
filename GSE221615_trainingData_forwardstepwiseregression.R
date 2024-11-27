library(lmtest)
library(ggplot2)
library(dplyr)
library(glmnet)

setwd("E:/Atherosclerosis_predictionModel_2024")
total_trainingdata_geneExpression <- read.csv("E:/Atherosclerosis_predictionModel_2024/GSE221615_trainingdata_geneExpression_w_clinicalFactors.csv", row.names = 1)
total_trainingdata_meta <- read.csv("E:/Atherosclerosis_predictionModel_2024/GSE221615_trainingdata_metadata.csv")

# null_test <- lm(total_trainingdata_geneExpression~1, data = total_trainingdata_geneExpression)
# lm(pesa_score_responseVariable~1, data = total_trainingdata_geneExpression)



###forward stepwise regression --> used to pick significant features
#the empty model would be a dataset with no coefficients (just pesa score)
##the empty model can be either a numeric or character vector
## data must contain BOTH dependent (pesa score) AND independent (gene ID, sex, age) 

pesa_score_responseVariable <- total_trainingdata_meta$pesa_score

total_trainingdata_geneExpression$pesa_score <- pesa_score_responseVariable
total_trainingdata_geneExpression <- total_trainingdata_geneExpression %>% select(pesa_score, everything())

##### not working -- stepwise regression#######
empty_forward_model <- lm(pesa_score ~ 1, data = total_trainingdata_geneExpression) 
all(is.na(total_trainingdata_geneExpression$pesa_score))

all(is.nan(total_trainingdata_geneExpression$pesa_score))
