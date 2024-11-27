library(glmnet)
library(dplyr)


setwd("E:/Atherosclerosis_predictionModel_2024/GSE221615_trainingData")

dataset <- read.csv("matrix_fpkm_GSE221615.csv", stringsAsFactors = F)
rownames(dataset) <- dataset$X
dataset <- dataset %>% select(-X)

0.8 * (388 - 27) ## 288.8
0.8 * 27 ## 21.6

metadata <- read.csv("metadata_GSE221615.csv", stringsAsFactors = F)

full_female_dataset <- metadata[metadata$gender == "female", ]
full_male_dataset <- metadata[metadata$gender == "male", ]

training_data_femaleSamples <- sample(full_female_dataset$Library.Name, size = 22)
training_data_maleSamples <- sample(full_male_dataset$Library.Name, size = 289)

training_data_femaleSamples <- as.data.frame(training_data_femaleSamples)
training_data_maleSamples <- as.data.frame(training_data_maleSamples)
metadata[metadata$Library.Name %in% training_data_femaleSamples$training_data_femaleSamples, ]

training_data_femaleSamples <- metadata[metadata$Library.Name %in% training_data_femaleSamples$training_data_femaleSamples, ]
training_data_maleSamples <-  metadata[metadata$Library.Name %in% training_data_maleSamples$training_data_maleSamples, ]

training_data_femaleSamples$Library.Name
trained_sampleID <- c(training_data_femaleSamples$Library.Name, training_data_maleSamples$Library.Name)

metadata[metadata$Library.Name %in% trained_sampleID, ]

training_dataset <- metadata[metadata$Library.Name %in% trained_sampleID, ]

testing_data_female <- full_female_dataset[!(full_female_dataset$Library.Name %in% training_data_femaleSamples$Library.Name), ]

testing_data <- metadata[!(metadata$Library.Name %in% training_dataset$Library.Name), ]
write.csv(testing_data, "metadata_testing_data_GSE221615.csv")
write.csv(training_dataset, "metadata_training_data_GSE221615.csv")

geneExpression_training_data <- dataset[, (colnames(dataset) %in% training_dataset$Library.Name) == "TRUE"]
colnames(dataset) %in% training_dataset$Library.Name
dataset[(colnames(dataset) %in% training_dataset$Library.Name) == "TRUE", ]
colnames(geneExpression_training_data) %in% training_dataset$Library.Name

testing_data <- setdiff(metadata$Library.Name, training_dataset$Library.Name)
testing_data <- as.data.frame(testing_data)
testing_data <- metadata[metadata$Library.Name %in% testing_data$testing_data, ]
geneExpression_testing_data <- dataset[, (colnames(dataset) %in% testing_data$Library.Name) == "TRUE"]

training_dataset <- training_dataset[order(match(training_dataset$Library.Name, colnames(geneExpression_training_data))), ]

disease_severity_responseVariable <- training_dataset$pesa_score
training_geneExpression_matrix <- t(geneExpression_training_data)
str(training_geneExpression_matrix)

cv.glmnet(training_geneExpression_matrix, disease_severity_responseVariable)
training_dataset <- training_dataset[training_dataset$Library.Name %in% rownames(training_geneExpression_matrix),]
disease_severity_responseVariable <- training_dataset$pesa_score
disesae_severity_trainingData_cv_Fit <- cv.glmnet(training_geneExpression_matrix, disease_severity_responseVariable, family = "multinomial")
plot(disesae_severity_trainingData_cv_Fit)
title("Multinomial Family Disease Severity on Training Dataset")

predict(disesae_severity_trainingData_cv_Fit, s = 0.060866, type = "nonzero")
x <- glmnet(training_geneExpression_matrix, disease_severity_responseVariable, lambda = 0.06680, family = "multinomial")
y <- as.data.frame(x)
y <- as.data.frame(x$df)
y$SE <- x$dim

library(coefplot)
extract.coef(x)




##sex age covarients 


copy_geneExpression_training_data <- geneExpression_training_data
training_dataset$numeric_sex <- training_dataset$gender
training_dataset$numeric_sex <- ifelse(training_dataset$numeric_sex=="male", 1, 0)
####male == 1, female = 0
copy_geneExpression_training_data <- rbind(copy_geneExpression_training_data, training_dataset$numeric_sex)
rownames(copy_geneExpression_training_data)
row.names(copy_geneExpression_training_data)[39377] <- "numeric_sex"
rownames(copy_geneExpression_training_data)
matrix_training_data_withSexNumeric <- as.matrix(copy_geneExpression_training_data)

matrix_training_data_withSexNumeric <- t(copy_geneExpression_training_data)
str(matrix_training_data_withSexNumeric)

cv.glmnet(matrix_training_data_withSexNumeric, disease_severity_responseVariable, family = "multinomial")
disesae_severity_sexcovarient_trainingData_cv_Fit <- cv.glmnet(matrix_training_data_withSexNumeric, disease_severity_responseVariable, family = "multinomial")
plot(disesae_severity_sexcovarient_trainingData_cv_Fit)
nonzero_coefs_sexcovariate <- predict(disesae_severity_sexcovarient_trainingData_cv_Fit, s = 0.07187, type = "nonzero")

focal_coefs_sexcovariate <- nonzero_coefs_sexcovariate$Focal
generalized_coefs_sexcovariate <- nonzero_coefs_sexcovariate$Generalized
control_coefs_sexcovariate <- nonzero_coefs_sexcovariate$No
intermediate_coef_sex_covariate <- nonzero_coefs_sexcovariate$Intermediate
annotated_generalized_coefs_sexcovariate <- Human.GRCh38.p13.annot.tsv[Human.GRCh38.p13.annot.tsv$GeneID %in% generalized_coefs_sexcovariate$X1, ]
annotated_no_coefs_sexcovariate <- Human.GRCh38.p13.annot.tsv[Human.GRCh38.p13.annot.tsv$GeneID %in% control_coefs_sexcovariate$X1, ]
annotated_intermediate_coefs_sexcovariate <- Human.GRCh38.p13.annot.tsv[Human.GRCh38.p13.annot.tsv$GeneID %in% intermediate_coef_sex_covariate$X1, ]


myCoefs <- coef(disesae_severity_sexcovarient_trainingData_cv_Fit, s=0.06998);
myCoefs[which(myCoefs$Generalized != 0 ) ]    

library(glmpath.cr)
