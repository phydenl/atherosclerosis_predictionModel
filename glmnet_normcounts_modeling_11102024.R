#####working the the submitter supplied normalized counts 
#####FPKM 
##likelihoood ratio test + relaxed lasso 
###extract nonzero genes 
###build regression model with selected genes + sex + age 
library(glmnet)
library(dplyr)
library(readr)

setwd("E:/Atherosclerosis_predictionModel_2024")

normcounts <- read_tsv("C:/Users/phyde/Downloads/GSE221615_normalizedCounts.tsv")
metadata <- read.csv("E:/Atherosclerosis_predictionModel_2024/GSE221615_trainingData/metadata_GSE221615.csv")
suaz <- read.csv("C:/Users/phyde/Downloads/SUAZ_to_GSM_sampleID - Sheet1.csv", header = F)


set.seed(42) ###did this before sample aka randomly selecting IDs from dataset 

colnames(normcounts)[3:length(colnames(normcounts))] <- suaz$V1

#####seperating dataset into male and female samples 
female_dataset <- normcounts[, 1:2]
male_dataset <- normcounts[,1:2]


female_metadata <- metadata[metadata$gender == "female", ]
male_metadata <- metadata[metadata$gender == "male", ]

colnames(normcounts)

female_dataset <- normcounts[, colnames(normcounts) %in% female_metadata$Library.Name ]
female_dataset$ens_gene <- normcounts$ens_gene
female_dataset <- female_dataset %>% select(ens_gene, everything())
male_dataset <- normcounts[, colnames(normcounts) %in% male_metadata$Library.Name]
male_dataset$ens_gene <- normcounts$ens_gene
male_dataset <- male_dataset %>% select(ens_gene, everything())

all(male_metadata$Library.Name %in% colnames(male_dataset)) ## true
all(female_metadata$Library.Name %in% colnames(female_dataset)) ## true

# write.csv(female_dataset, "NormCounts_GSE221615_femaleGeneExpression_total.csv")
# write.csv(male_dataset, "NormCounts_GSE221615_maleGeneExpression_total.csv")

####selecting a random 80% to be training data 
0.8*27 ##21.6
0.8 * 364 ##291.2 
female_training_metadata <- sample(female_metadata$Library.Name, size = 22)
female_training_metadata <- as.data.frame(female_training_metadata)
colnames(female_training_metadata) <- "Library.Name"
female_training_metadata <- right_join(female_metadata, female_training_metadata, by = "Library.Name")

male_training_metadata <- as.data.frame(sample(male_metadata$Library.Name, size = 291))
colnames(male_training_metadata) <- "Library.Name"
male_training_metadata <- right_join(male_metadata, male_training_metadata, by = "Library.Name")


male_geneExpression_trainingdata <- normcounts[, colnames(normcounts) %in% male_training_metadata$Library.Name]
female_geneExpression_trainingdata <- normcounts[, colnames(normcounts) %in% female_training_metadata$Library.Name]

# write.csv(female_training_metadata, "normcounts_GSE221615_femaleTrainingData_metadata.csv")
# write.csv(male_training_metadata, "normcounts_GSE221615_maleTrainingData_metadata.csv")
# write.csv(male_geneExpression_trainingdata, "normcounts_GSE221615_maleTrainingData_geneExpression.csv")
# write.csv(female_geneExpression_trainingdata, "normcounts_GSE221615_femaleTrainingData_geneExpression.csv")

## making ens_gene the rownames for female and male gene expression training dataset
##rownames(female_geneExpression_trainingdata) <- normcounts$ens_gene
##rownames(male_geneExpression_trainingdata) <- normcounts$ens_gene

## making genes as columns and sample IDs as rownames using t() function 
# matrix_female_geneExpression_trainingdata <- t(female_geneExpression_trainingdata)
# matrix_female_geneExpression_trainingdata<- as.data.frame(matrix_female_geneExpression_trainingdata)
# matrix_female_geneExpression_trainingdata$numeric_sex <- 0
# 
# matrix_male_geneExpression_trainingdata <- t(male_geneExpression_trainingdata)
# matrix_male_geneExpression_trainingdata <- as.data.frame(matrix_male_geneExpression_trainingdata)
# matrix_male_geneExpression_trainingdata$numeric_sex <- 1 
# 
# 
# female_geneExpression_trainingdata_numericsex <- as.data.frame(t(matrix_female_geneExpression_trainingdata))
# male_geneExpression_trainingdata_numericsex <- as.data.frame(t(matrix_male_geneExpression_trainingdata))
# 
# 
# test <- merge(female_geneExpression_trainingdata_numericsex, male_geneExpression_trainingdata_numericsex)



total_trainingdata <- metadata[metadata$Library.Name %in% female_training_metadata$Library.Name | 
                                 metadata$Library.Name %in% male_training_metadata$Library.Name, ]

total_trainingdata_geneExpression <- normcounts[, colnames(normcounts) %in% total_trainingdata$Library.Name, ]
total_trainingdata$numeric_sex <- total_trainingdata$gender
total_trainingdata$numeric_sex <- ifelse(total_trainingdata$numeric_sex=="male", 1, 0)
rownames(total_trainingdata_geneExpression) <- normcounts$ens_gene
total_trainingdata_geneExpression <- as.data.frame(t(total_trainingdata_geneExpression))

total_trainingdata <- total_trainingdata[order(match(total_trainingdata$Library.Name, rownames(total_trainingdata_geneExpression))), ]
total_trainingdata_geneExpression$numeric_sex <- NA
total_trainingdata_geneExpression$numeric_sex <- total_trainingdata$numeric_sex
total_trainingdata_geneExpression$age <- total_trainingdata$Age


#### getting lambda multinomial regression from cv.glmnet
##genes on top, sample ID on the side 
set.seed(42)
pesa_score_responseVariable <- total_trainingdata$pesa_score
matrix_total_trainingdata_geneExpression <- as.matrix(total_trainingdata_geneExpression)

cv_fit_trainingdata_w_clinicalfactors <- cv.glmnet(matrix_total_trainingdata_geneExpression, pesa_score_responseVariable, family = "multinomial", nfolds=5)

cv_fit_trainingdata_w_clinicalfactors$lambda.min
predict(cv_fit_trainingdata_w_clinicalfactors, s =cv_fit_trainingdata_w_clinicalfactors$lambda.min, type = "nonzero")
withsex_min_lambda <- cv_fit_trainingdata_w_clinicalfactors$lambda.min
non_zero_coefs <- coef(cv_fit_trainingdata_w_clinicalfactors, s =withsex_min_lambda )
#non_zero_genes <- colnames(matrix_total_trainingdata_geneExpression)[which(non_zero_coefs != 0)]
non_zero_coefs_generalized <- non_zero_coefs$Generalized
dimnames(non_zero_coefs_generalized)[[1]]
non_zero_coefs_generalized@x

plot(cv_fit_trainingdata_w_clinicalfactors)

#row_indicies <- row(non_zero_coefs_generalized)[non_zero_coefs_generalized != 0]
#tmp_coeffs <- coef(cv.glmnet.fit, s = "lambda.min")
non_zero_generalized_trainingdata <- data.frame(name = non_zero_coefs_generalized@Dimnames[[1]][non_zero_coefs_generalized@i + 1], coefficient = non_zero_coefs_generalized@x)
non_zero_intermediate_trainingdata<- data.frame(name = non_zero_coefs$Intermediate@Dimnames[[1]][non_zero_coefs$Intermediate@i + 1], coefficient = non_zero_coefs$Intermediate@x)
non_zero_control_trainingdata <- data.frame(name = non_zero_coefs$No@Dimnames[[1]][non_zero_coefs$No@i + 1], coefficient = non_zero_coefs$No@x)
non_zero_focal_trainingdata <- data.frame(name = non_zero_coefs$Focal@Dimnames[[1]][non_zero_coefs$Focal@i + 1], coefficient = non_zero_coefs$Focal@x)

# write.csv(non_zero_generalized_trainingdata, "GSE221615_trainingdata_nonzerocoef_generalized_withsex.csv")
# write.csv(non_zero_intermediate_trainingdata, "GSE221615_trainingdata_nonzerocoef_intermediate_withsex.csv")
# write.csv(non_zero_focal_trainingdata, "GSE221615_trainingdata_nonzerocoef_focal_withsex.csv")
# write.csv(non_zero_control_trainingdata, "GSE221615_trainingdata_nonzerocoef_control_withsex.csv")


total_trainingdata_geneexpression_no_clinical <- total_trainingdata_geneExpression %>% select(-"numeric_sex")
set.seed(42)
total_trainingdata_geneexpression_no_clinical_matrix <- as.matrix(total_trainingdata_geneexpression_no_clinical)
total_trainingdata_geneexpression_no_clinical_cvfit <- cv.glmnet(total_trainingdata_geneexpression_no_clinical_matrix, pesa_score_responseVariable, family = "multinomial", nfolds = 5)

min_lambda <- total_trainingdata_geneexpression_no_clinical_cvfit$lambda.min
plot(total_trainingdata_geneexpression_no_clinical_cvfit)
fit_nofactors <- glmnet(total_trainingdata_geneexpression_no_clinical_matrix, pesa_score_responseVariable, lambda = min_lambda, family = "multinomial")

fit_nofactors_generalized <- data.frame(name = fit_nofactors$beta$Generalized@Dimnames[[1]][fit_nofactors$beta$Generalized@i + 1], coefficient = fit_nofactors$beta$Generalized@x)
fit_nofactors_intermediate<- data.frame(name = fit_nofactors$beta$Intermediate@Dimnames[[1]][fit_nofactors$beta$Intermediate@i + 1], coefficient = fit_nofactors$beta$Intermediate@x)
fit_nofactors_control <- data.frame(name = fit_nofactors$beta$No@Dimnames[[1]][fit_nofactors$beta$No@i + 1], coefficient = fit_nofactors$beta$No@x)
fit_nofactors_focal <- data.frame(name =  fit_nofactors$beta$Focal@Dimnames[[1]][fit_nofactors$beta$Focal@i + 1], coefficient = fit_nofactors$beta$Focal@x)

# write.csv(fit_nosex_control, "GSE221615_trainingdata_nonzerocoef_control_withoutsex.csv")
# write.csv(fit_nosex_generalized, "GSE221615_trainingdata_nonzerocoef_generalized_withoutsex.csv")


############relaxed lasso


#########################################
###     RELAXED LASSO NOT WORKING!    ##
#########################################
# 
# fit_w_clinicalfactors <- glmnet(matrix_total_trainingdata_geneExpression, pesa_score_responseVariable, relax = T, family = "multinomial")
# ##dimension error 
# dim(matrix_total_trainingdata_geneExpression)
# length(pesa_score_responseVariable)
# levels(pesa_score_responseVariable) #### pesa-score is treated as "NULL", use as.factor 
# pesa_score_responseVariable <- as.factor(pesa_score_responseVariable)
# levels(pesa_score_responseVariable)
# ##should be fixed
# fit_w_clinicalfactors <- glmnet(matrix_total_trainingdata_geneExpression, pesa_score_responseVariable, relax = T, family = "multinomial")
# fit_w_clinicalfactors <- glmnet(matrix_total_trainingdata_geneExpression, 
#                       pesa_score_responseVariable, 
#                       family = "multinomial", 
#                       maxit = 200000, 
#                       relax = TRUE)
# 
# str(pesa_score_responseVariable)
# fit_withsex <- glmnet(matrix_total_trainingdata_geneExpression, pesa_score_responseVariable, relax = T, family = "multinomial")
# dim(matrix_total_trainingdata_geneExpression)
# length(pesa_score_responseVariable)
# # Scale the data if needed
# scaled_matrix_total_trainingdata_geneExpression <- scale(matrix_total_trainingdata_geneExpression)
# fit_withsex <- glmnet(scaled_matrix_total_trainingdata_geneExpression, pesa_score_responseVariable, relax = T, family = "multinomial", maxit = 1000000)
# pesa_score_response <- as.factor(pesa_score_responseVariable)
# unique(levels(pesa_score_response))
# any(is.na(pesa_score_response))
# any(is.na(pesa_score_responseVariable))
# any(is.na(matrix_total_trainingdata_geneExpression))
# 
# fit_withsex_relaxed <- glmnet(matrix_total_trainingdata_geneExpression, 
#                       pesa_score_responseVariable, 
#                       family = "multinomial", 
#                       relax = TRUE, 
#                       maxit = 200000, 
#                       lambda = 0.0001)
# 
# fit_withsexrelaxed_generalized <- data.frame(name = fit_withsex_relaxed$beta$Generalized@Dimnames[[1]][fit_withsex_relaxed$beta$Generalized@i + 1], 
#                                              coefficient = fit_withsex_relaxed$beta$Generalized@x)
# fit_withsexrelaxed_intermediate<- data.frame(name = fit_withsex_relaxed$beta$Intermediate@Dimnames[[1]][fit_withsex_relaxed$beta$Intermediate@i + 1], coefficient = fit_withsex_relaxed$beta$Intermediate@x)
# fit_withsexrelaxed_control <- data.frame(name = fit_withsex_relaxed$beta$No@Dimnames[[1]][fit_withsex_relaxed$beta$No@i + 1], coefficient = fit_withsex_relaxed$beta$No@x)
# fit_withsexrelaxed_focal <- data.frame(name =  fit_withsex_relaxed$beta$Focal@Dimnames[[1]][fit_withsex_relaxed$beta$Focal@i + 1], coefficient = fit_withsex_relaxed$beta$Focal@x)
# 
# 
# fit_withoutsex_relaxed <- cv.glmnet(total_trainingdata_geneexpression_nosex, 
#                                  pesa_score_responseVariable, 
#                                  family = "multinomial", 
#                                  relax = TRUE, 
#                                  gamma = 0,
#                                  maxp=150, 
#                                  path = T)
# plot(fit_withoutsex_relaxed)
# 
# coef_matrix <- coef(fit_withoutsex_relaxed)
# 
# # Extract non-zero coefficients
# non_zero_coefs <- coef_matrix[, length(fit_withoutsex_relaxed$lambda)]  # for the final lambda
# non_zero_coefs <- non_zero_coefs[non_zero_coefs != 0]  # filter out zeros
# fit_withoutsex_relaxed$
# fit_nosex_generalized <- data.frame(name = fit_withoutsex_relaxed$beta$Generalized@Dimnames[[1]][fit_withoutsex_relaxed$beta$Generalized@i + 1], 
#                                     coefficient = fit_withoutsex_relaxed$beta$Generalized@x)
# fit_nosex_intermediate<- data.frame(name = fit_withoutsex_relaxed$beta$Intermediate@Dimnames[[1]][fit_withoutsex_relaxed$beta$Intermediate@i + 1], 
#                                     coefficient = fit_withoutsex_relaxed$beta$Intermediate@x)
# fit_nosex_control <- data.frame(name = fit_withoutsex_relaxed$beta$No@Dimnames[[1]][fit_withoutsex_relaxed$beta$No@i + 1], 
#                                 coefficient = fit_withoutsex_relaxed$beta$No@x)
# fit_nosex_focal <- data.frame(name =  fit_withoutsex_relaxed$beta$Focal@Dimnames[[1]][fit_withoutsex_relaxed$beta$Focal@i + 1], 
#                               coefficient = fit_withoutsex_relaxed$beta$Focal@x)


# write.csv(fit_withsexrelaxed_generalized, "GSE221615_trainingdata_nonzerocoef_generalized_withsex_relaxedlasso.csv")
# write.csv(fit_withsexrelaxed_intermediate, "GSE221615_trainingdata_nonzerocoef_intermediate_withsex_relaxedlasso.csv")
# write.csv(fit_withsexrelaxed_control, "GSE221615_trainingdata_nonzerocoef_focal_withsex_relaxedlasso.csv")
# write.csv(fit_withsexrelaxed_focal, "GSE221615_trainingdata_nonzerocoef_control_withsex_relaxedlasso.csv")
# 
# 
# 
# write.csv(fit_nosex_generalized, "GSE221615_trainingdata_nonzerocoef_generalized_nosex_relaxedlasso.csv")
# write.csv(fit_nosex_intermediate, "GSE221615_trainingdata_nonzerocoef_intermediate_nosex_relaxedlasso.csv")
# write.csv(fit_nosex_control, "GSE221615_trainingdata_nonzerocoef_focal_nosex_relaxedlasso.csv")
# write.csv(fit_nosex_focal, "GSE221615_trainingdata_nonzerocoef_control_nosex_relaxedlasso.csv")
# 
# write.csv(total_trainingdata_geneExpression, "GSE221615_trainingdata_for_stepwiseRegression.csv")
# write.csv(total_trainingdata, "GSE221615_trainingdata_metadata.csv")
