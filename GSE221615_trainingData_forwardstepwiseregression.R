#install.packages("familiar", dependencies = TRUE)


library(lmtest)
library(ggplot2)
library(dplyr)
library(glmnet)
library(familiar)

# GLOBAL VARIABLES
# set this to the directory that contains the scripts and data
working_dir = "C:\\Users\\derek\\OneDrive - University of Connecticut\\software\\atherosclerosis_predictionModel\\" 


setwd(working_dir)
total_trainingdata_geneExpression <- read.csv("GSE221615_trainingdata_geneExpression_w_clinicalFactors.csv", row.names = 1)
total_trainingdata_meta <- read.csv("GSE221615_trainingdata_metadata.csv")

# null_test <- lm(total_trainingdata_geneExpression~1, data = total_trainingdata_geneExpression)
# lm(pesa_score_responseVariable~1, data = total_trainingdata_geneExpression)



###forward stepwise regression --> used to pick significant features
#the empty model would be a dataset with no coefficients (just pesa score)
##the empty model can be either a numeric or character vector
## data must contain BOTH dependent (pesa score) AND independent (gene ID, sex, age) 

pesa_score_responseVariable <- total_trainingdata_meta$pesa_score

total_trainingdata_geneExpression$pesa_score <- as.factor(pesa_score_responseVariable)
total_trainingdata_geneExpression <- total_trainingdata_geneExpression %>% select(pesa_score, everything())

##### not working -- stepwise regression#######
empty_forward_model <- glmnet(1, pesa_score, family = "multinomial", data = total_trainingdata_geneExpression) 

cv_fit_trainingdata_w_clinicalfactors <- cv.glmnet(matrix_total_trainingdata_geneExpression, pesa_score_responseVariable, family = "multinomial", nfolds=5)

familiar::summon_familiar(
  data = total_trainingdata_geneExpression,
  experiment_dir = file.path(working_dir, "output_of_familiar"),
  outcome_type = "multinomial",
  outcome_column = "pesa_score",
  experimental_design = "fs+mb",
  cluster_method = "none",
  fs_method = "lasso_multinomial",
  learner = "glm",
  parallel = FALSE)