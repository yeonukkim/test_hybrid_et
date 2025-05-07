rm(list = ls())

# load libraries
library(tidyverse)
library(caret)
library(parallel)
library(doParallel)

# load dataset
df <- read.csv("FLX_RS_merged_filtered2.csv") 
head(df)

# predictors
predictors <- c("Fpar","soil_moisture","VegWC","TA_F","CO2_F_MDS","SW_IN_F_MDS","ga_H","rha")
predictors2 <- c("Fpar","soil_moisture","VegWC","TA_F","CO2_F_MDS","SW_IN_F_MDS","ga_H","rha", "AE")

# 75% of data used for model tuning/validation
set.seed(123)
index <- createDataPartition(df[,1], p=0.75, list=F)
saveRDS(index, "data_partition_index.rds")
# index <- readRDS("data_partition_index.rds")

# target variable list
target_list <- c("LE_CORR", 
								 "PT_stress_cor", 
								 "Penman_stress_cor", 
								 "PM_stress_cor", 
								 "Frh_cor", 
								 "drh_cor", 
								 "Gs_log_cor")


##############################################
#### rf run for mtry tuning at ntree=100
##############################################

#Add parallel processing for the fast processing if you want to
cl <- makePSOCKcluster(5)
registerDoParallel(cl)

for(i in c(1:7)){
	# variable we need for predict 
	target <- target_list[i]
	print(target)
	
	if(i <= 1){
		ML.df <- df %>% select(c(target,predictors2))
		tgrid <- data.frame(mtry = c(1:9))
	} else{
		ML.df <- df %>% select(c(target,predictors))
		tgrid <- data.frame(mtry = c(1:8))
		}
	
	train_set <- ML.df[index,]
	test_set <- ML.df[-index,]
	
	############### Random forest run
	#### random forest model with mtry tuning
	
	set.seed(42)
	RF <- train(formula(paste(target, "~.")),
							data = train_set,
							method = "rf",
							preProcess = c("medianImpute"),                
							trControl=trainControl(method = "cv",           
																	 number = 3),
							tuneGrid = tgrid,
							na.action = na.pass,
							allowParallel=TRUE, 
							ntree=100, 
							importance = TRUE)
	
	saveRDS(RF, paste0("./ntree_models/RF_",target,"_ntree_100.rds"))
}

stopCluster(cl)

### mtry check
RF_LE_CORR <- readRDS("RF_LE_CORR.rds")
RF_PT_stress_cor <- readRDS("RF_PT_stress_cor.rds")
RF_Penman_stress_cor <- readRDS("RF_Penman_stress_cor.rds")
RF_PM_stress_cor <- readRDS("RF_PM_stress_cor.rds")
RF_Gs_log_cor <- readRDS("RF_Gs_log_cor.rds")
RF_Frh_cor <- readRDS("RF_Frh_cor.rds")
RF_drh_cor <- readRDS("RF_drh_cor.rds")

RF_LE_CORR$bestTune
RF_PT_stress_cor$bestTune
RF_Penman_stress_cor$bestTune
RF_PM_stress_cor$bestTune
RF_Gs_log_cor$bestTune
RF_Frh_cor$bestTune
RF_drh_cor$bestTune

##############################################
#### rf run with best mtry at ntree=500
##############################################
# best mtry when ntree = 100
mtry <- c(6,6,3,3,3,6,4)

#Add parallel processing for the fast processing if you want to
cl <- makePSOCKcluster(5)
registerDoParallel(cl)

for(i in c(1:7)){
  # variable we need for predict 
  target <- target_list[i]
  print(target)
  
  if(i == 1){
    ML.df <- df %>% select(c(target,predictors2))
  } else{
    ML.df <- df %>% select(c(target,predictors))
  }
  
  train_set <- ML.df[index,]
  test_set <- ML.df[-index,]
  
  ############### Random forest run
  set.seed(42)
  RF <- train(formula(paste(target, "~.")),
              data = train_set,
              method = "rf",
              preProcess = c("medianImpute"),                
              trControl=trainControl(method = "none"),
              tuneGrid = data.frame(mtry = mtry[i]),
              na.action = na.pass,
              allowParallel=TRUE, 
              ntree=500, 
              importance = TRUE)
  
  saveRDS(RF, paste0("./ntree_models/RF_",target,"_ntree_500.rds"))
    
}

stopCluster(cl)

