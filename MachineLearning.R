# New performance sampling, alle experiments in one file
require(caret)
require(randomForest)
require(pROC)
require(doMC)
require(rjson)

nExperiments = 100
nCores = 20
registerDoMC(cores = nCores)


#TODO: Identify which feature set you are using, for the experiments with the provenance and the predicate features.
## Guney reference set
#PredColumns = c(2:142)
#ProvColumns = c(143:213)

# EMC reference set
PredColumns = c(2:149)
ProvColumns = c(150:227)

# General code
## Read in the data
f = read.csv2("Feature set.csv", stringsAsFactors = F)
input = f

reference_set = read.csv2("Reference set.csv", stringsAsFactors = F)
colnames(reference_set) = c("Drug", "DrugTargets")

## Calculate the number of drug targets per drug, and add this information as a feature
reference_set$nDT = sapply(reference_set$DrugTargets, function(x){length(unique(fromJSON(x)))})
input = merge(input, unique(reference_set[,c("Drug", "nDT")]), by = "Drug")

## Separate the dataset into the features and the none-features
removalColumns = c("Drug", "Disease", "diseaseProtCount", "nDT")
dd_data = input[,which(colnames(input) %in% removalColumns)]
input = input[,-which(colnames(input) %in% removalColumns)]

## Function for sampling of the instances
sample_count = length(which(input$goldstandard == "VALID"))
createBalancedSet = function(df, n, pos_class = "VALID"){
  balanced = rbind(df[sample(row.names(df[df$goldstandard == pos_class, ]), n),], 
                   df[sample(row.names(df[!df$goldstandard == pos_class, ]), n),])
  return(balanced)}

# Create the data frame which is used as output
out = data.frame( auc_complete = numeric(nExperiments),
                  auc_oneTarget = numeric(nExperiments),
                  auc_twoTarget = numeric(nExperiments),
                  auc_threeTarget = numeric(nExperiments),
                  auc_indirect_only = numeric(nExperiments),
                  auc_preds_only = numeric(nExperiments),
                  auc_prov_only = numeric(nExperiments),
                  auc_noFeatures = numeric(nExperiments))

## Calculate the performances
for(i in 1:nrow(out)){
  print(i)
  registerDoMC(cores = nCores)
  input_sample = createBalancedSet(input, sample_count)
  dd_sample = input_sample[,which(colnames(input_sample) %in% removalColumns)]
  dd_sample$rowIndex = 1:nrow(input_sample)
  input_sample = input_sample[,-which(colnames(input_sample) %in% removalColumns)]
  
  # Complete dataset
  model = train(goldstandard ~., data = input_sample, method = "rf",
                trControl=trainControl(method="cv",number=10, classProbs = TRUE, 
                                       savePredictions = TRUE, summaryFunction = twoClassSummary),
                prox=TRUE, allowParallel=TRUE, metric = "ROC")
  out$auc_complete[i] = model$results$ROC[model$results$mtry == model$bestTune$mtry]
  
  # See where the drugs with 1, 2, or 3 or more drugtargets are classified
  cv_results = model$pred[model$pred$mtry == model$bestTune$mtry, ]
  cv_results = merge(cv_results[,c("obs", "VALID", "rowIndex")], dd_sample[,c("nDT", "rowIndex")], by = "rowIndex")
  out$auc_oneTarget[i] = as.numeric(roc(cv_results$obs[cv_results$nDT == 1], cv_results$VALID[cv_results$nDT == 1])$auc)
  out$auc_twoTarget[i] = as.numeric(roc(cv_results$obs[cv_results$nDT == 2], cv_results$VALID[cv_results$nDT == 2])$auc)
  out$auc_threeTarget[i] = as.numeric(roc(cv_results$obs[cv_results$nDT > 2], cv_results$VALID[cv_results$nDT > 2])$auc)
  
  # Remove the direct relationships and measure performance
  input_nobias = input_sample[,-grep("direct_", colnames(input_sample))]
  input_nobias = input_nobias[,-grep("selfref_", colnames(input_nobias))]
  input_nobias$overlap = NULL
  
  model = train(goldstandard ~., data = input_nobias, method = "rf",
                trControl=trainControl(method="cv",number=10, classProbs = TRUE, 
                                       savePredictions = TRUE, summaryFunction = twoClassSummary),
                prox=TRUE, allowParallel=TRUE, metric = "ROC")
  out$auc_indirect_only[i] = model$results$ROC[model$results$mtry == model$bestTune$mtry]
  
  # Remove the provenances and train
  input_noProv = input_sample[,-ProvColumns]
  model = train(goldstandard ~., data = input_noProv, method = "rf",
                  trControl=trainControl(method="cv",number=10, classProbs = TRUE, 
                                         savePredictions = TRUE, summaryFunction = twoClassSummary),
                  prox=TRUE, allowParallel=TRUE, metric = "ROC")
  out$auc_preds_only[i] = model$results$ROC[model$results$mtry == model$bestTune$mtry]  
    
  # Remove the predicates and train
  input_noPred = input_sample[,-PredColumns]
  
  model = train(goldstandard ~., data = input_noPred, method = "rf",
                trControl=trainControl(method="cv",number=10, classProbs = TRUE, 
                                       savePredictions = TRUE, summaryFunction = twoClassSummary),
                prox=TRUE, allowParallel=TRUE, metric = "ROC")
  out$auc_prov_only[i] = model$results$ROC[model$results$mtry == model$bestTune$mtry] 
  
  # Create only "Associated with" relationships and train
  direct_cols = grep("direct_", colnames(input_sample))
  indirect_cols = grep(c("step1_|step2_"), colnames(input_sample))
  
  input_sample$direct_associated_with = ifelse(rowSums(input_sample[,direct_cols]) > 0, TRUE, FALSE)
  input_sample$indirect_associated_with = ifelse(rowSums(input_sample[,indirect_cols]) > 0, TRUE, FALSE)
  
  input_sample = input_sample[,c("overlap", "direct_associated_with", "indirect_associated_with", "goldstandard")]

  model = train(goldstandard ~., data = input_sample, method = "rf",
                trControl=trainControl(method="cv",number=10, classProbs = TRUE, 
                                       savePredictions = TRUE, summaryFunction = twoClassSummary),
                prox=TRUE, allowParallel=TRUE, metric = "ROC")
  out$auc_noFeatures[i] = model$results$ROC[model$results$mtry == model$bestTune$mtry]
}

write.csv2(out, "Results crossvalidation experiments.csv", row.names = F)