#### run random forest to find the upper limit of the intermediate parameters
#### in some hybrid models. Select upper limit which yield the lowest rmse.
#### Note: hybrid 5 and 6 models do not require this analysis.

rm(list = ls())

# load libraries
library(tidyverse)
library(caret)
library(parallel)
library(doParallel)
library(Metrics)
library(bigleaf)

# load dataset
df <- read.csv("FLX_RS_merged_filtered2.csv") 
head(df)

# predictors
predictors <- c("Fpar","soil_moisture","VegWC","TA_F","CO2_F_MDS","SW_IN_F_MDS","ga_H","rha")

# only use training set in this preliminary run
set.seed(123)
index <- createDataPartition(df[,1], p=0.75, list=F)
df_train <- df[index,]

# define data frame
results <- data.frame(percentile = c(99,99.05,99.1,99.15,99.2,99.25,99.3,99.35,99.4,99.45,99.5,99.55,99.6,99.65,99.7,99.75,99.8,99.85,99.9,99.95),
                      PT_max = rep(NA,20), PT_rmse = rep(NA,20),
                      Penman_max = rep(NA,20), Penman_rmse = rep(NA,20),
                      PMKc_max = rep(NA,20), PMKc_rmse = rep(NA,20),
                      PMgs_max = rep(NA,20), PMgs_min = rep(NA,20), PMgs_rmse = rep(NA,20)
)

###########################
#### model run
###########################
df_test <- df[-index,]

#Add parallel processing
cl <- makePSOCKcluster(5)
registerDoParallel(cl)

### PT 
for(j in c(1:20)){
	i <- results$percentile[j]
	# variable we need for predict 
	target <- "PT_stress_cor_raw"
	
	ML.df <- df_train %>% select(c(target,predictors))
	
	set.seed(123)
	index2 <- createDataPartition(ML.df[,1], p=0.75, list=F) 
	
	train_set <- ML.df[index2,]
	valid_set <- df_train[-index2,]
	
	use <- df %>% filter(AE > 30) %>%
		  	summarise(CORR = quantile(LE_CORR/PET_PT,i/100,na.rm=T))
	
	train_set$PT_stress_cor_raw <- ifelse(train_set$PT_stress_cor_raw > use$CORR[1], use$CORR[1], train_set$PT_stress_cor_raw)
	
	############### Random forest run
	#### random forest model with default mtry
	tgrid <- data.frame(mtry = 3)
	
	set.seed(42)
	RF <- train(formula(paste(target, "~.")),
							data = train_set,
							method = "rf",
							preProcess = c("medianImpute"),               
							trControl=trainControl(method = "none"),
							tuneGrid = tgrid,
							na.action = na.pass,
							allowParallel=TRUE, 
							ntree=100, 
							importance = FALSE)
	
	valid_set$PT_stress_cor_RF <- predict(RF, valid_set, na.action = na.pass)
	valid_set$LE_CORR_model <- valid_set$PT_stress_cor_RF * valid_set$PET_PT
	
	print(i)
	print(use$CORR[1])
	print(rmse(valid_set$LE_CORR_model,valid_set$LE_CORR))
	
	results$PT_max[j] <- use$CORR[1]
	results$PT_rmse[j] <- rmse(valid_set$LE_CORR_model,valid_set$LE_CORR)
}


### Penman
for(j in c(1:20)){
	i <- results$percentile[j]
	# variable we need for predict 
	target <- "Penman_stress_cor_raw"
	
	ML.df <- df_train %>% select(c(target,predictors))
	
	set.seed(123)
	index2 <- createDataPartition(ML.df[,1], p=0.75, list=F) 
	
	train_set <- ML.df[index2,]
	valid_set <- df_train[-index2,]
	
	use <- df %>% filter(AE > 30) %>%
		summarise(CORR = quantile(LE_CORR/PET_Penman,i/100,na.rm=T))
	
	train_set$Penman_stress_cor_raw <- ifelse(train_set$Penman_stress_cor_raw > use$CORR[1], use$CORR[1], train_set$Penman_stress_cor_raw)
	
	############### Random forest run
	#### random forest model with default mtry
	tgrid <- data.frame(mtry = 3)
	
	set.seed(42)
	RF <- train(formula(paste(target, "~.")),
							data = train_set,
							method = "rf",
							preProcess = c("medianImpute"),                
							trControl=trainControl(method = "none"),
							tuneGrid = tgrid,
							na.action = na.pass,
							allowParallel=TRUE,
							ntree=100, 
							importance = FALSE)
	
	valid_set$Penman_stress_cor_RF <- predict(RF, valid_set, na.action = na.pass)
	valid_set$LE_CORR_model <- valid_set$Penman_stress_cor_RF * valid_set$PET_Penman
	
	print(i)
	print(use$CORR[1])
	print(rmse(valid_set$LE_CORR_model,valid_set$LE_CORR))
	
	results$Penman_max[j] <- use$CORR[1]
	results$Penman_rmse[j] <- rmse(valid_set$LE_CORR_model,valid_set$LE_CORR)
}


### PM-Kc 
for(j in c(1:20)){
	i <- results$percentile[j]
	# variable we need for predict 
	target <- "PM_stress_cor_raw"
	
	ML.df <- df_train %>% select(c(target,predictors))
	
	set.seed(123)
	index2 <- createDataPartition(ML.df[,1], p=0.75, list=F) 
	
	train_set <- ML.df[index2,]
	valid_set <- df_train[-index2,]
	
	use <- df %>% filter(AE > 30) %>%
		summarise(CORR = quantile(LE_CORR/PET_PM,i/100,na.rm=T))
	
	train_set$PM_stress_cor_raw <- ifelse(train_set$PM_stress_cor_raw > use$CORR[1], use$CORR[1], train_set$PM_stress_cor_raw)
	
	############### Random forest run
	#### random forest model with default mtry
	tgrid <- data.frame(mtry = 3)
	
	set.seed(42)
	RF <- train(formula(paste(target, "~.")),
							data = train_set,
							method = "rf",
							preProcess = c("medianImpute"),                
							trControl=trainControl(method = "none"),
							tuneGrid = tgrid,
							na.action = na.pass,
							allowParallel=TRUE, 
							ntree=100, 
							importance = FALSE)
	
	valid_set$PM_stress_cor_RF <- predict(RF, valid_set, na.action = na.pass)
	valid_set$LE_CORR_model <- valid_set$PM_stress_cor_RF * valid_set$PET_PM
	
	print(i)
	print(use$CORR[1])
	print(rmse(valid_set$LE_CORR_model,valid_set$LE_CORR))
	
	results$PMKc_max[j] <- use$CORR[1]
	results$PMKc_rmse[j] <- rmse(valid_set$LE_CORR_model,valid_set$LE_CORR)
}



### PM-gs 
Gs_cor2 <- df$Gs_cor_raw

for(j in c(1:20)){
	i <- results$percentile[j]
	# variable we need for predict 
	target <- "Gs_log_cor"
	
	df_train <- df_train %>%
		mutate(Gs_cor = ifelse(!is.na(Gs_cor_raw), Gs_cor_raw,
													 ifelse(LE_CORR <= 0, quantile(Gs_cor2,1-i/100,na.rm=T),
													 			 ifelse(is.infinite(Gs_cor_raw),quantile(Gs_cor2,1-i/100,na.rm=T),quantile(Gs_cor2,i/100,na.rm=T)))),
					 Gs_log_cor = log(Gs_cor)
		)
	
	ML.df <- df_train %>% select(c(target,predictors))
	
	set.seed(123)
	index2 <- createDataPartition(ML.df[,1], p=0.75, list=F) 
	
	train_set <- ML.df[index2,]
	valid_set <- df_train[-index2,]
	
	
	############### Random forest run
	#### random forest model with default mtry
	tgrid <- data.frame(mtry = 3)
	
	set.seed(42)
	RF <- train(formula(paste(target, "~.")),
							data = train_set,
							method = "rf",
							preProcess = c("medianImpute"),                
							trControl=trainControl(method = "none"),
							tuneGrid = tgrid,
							na.action = na.pass,
							allowParallel=TRUE,
							ntree=100, 
							importance = FALSE)
	
	valid_set$Gs_log_cor_RF <- predict(RF, valid_set, na.action = na.pass)
	
	cp <- bigleaf.constants()$cp
	valid_set$density <- air.density(valid_set$TA_F,valid_set$PA_F)
	valid_set$esat_air <- Esat.slope(valid_set$TA_F)$Esat
	valid_set$Delta <- Esat.slope(valid_set$TA_F)$Delta
	valid_set$gamma <- psychrometric.constant(valid_set$TA_F,valid_set$PA_F)
	valid_set$Lv <- latent.heat.vaporization(valid_set$TA_F)
	
	valid_set$LE_CORR_model <- valid_set %>%
		mutate(LE_CORR_hybrid_PMgs = (Delta * (NETRAD - G_F_MDS) + density * cp * VPD_F / 10 / ra_H) / 
					 	(Delta + gamma + gamma /exp(Gs_log_cor_RF) / ra_H)) %>%
		pull(LE_CORR_hybrid_PMgs)
	
	print(i)
	print(quantile(Gs_cor2,1-i/100,na.rm=T))
	print(quantile(Gs_cor2,i/100,na.rm=T))
	print(rmse(valid_set$LE_CORR_model,valid_set$LE_CORR))
	
	results$PMgs_min[j] <- quantile(Gs_cor2,1-i/100,na.rm=T)
	results$PMgs_max[j] <- quantile(Gs_cor2,i/100,na.rm=T)
	results$PMgs_rmse[j] <- rmse(valid_set$LE_CORR_model,valid_set$LE_CORR)
	
}

stopCluster(cl)


write.csv(results,"max_decision.csv")

########################################################################
################ upper limit of the target parameter with lowest rmse 
########################################################################

results %>%
	ggplot(aes(percentile,PT_rmse)) +
	geom_point()

results %>%
	ggplot(aes(percentile,Penman_rmse)) +
	geom_point()

results %>%
	ggplot(aes(percentile,PMKc_rmse)) +
	geom_point()

results %>%
	ggplot(aes(percentile,PMgs_rmse)) +
	geom_point()

results %>%
	ggplot(aes(percentile,PMgs_max)) +
	geom_point()


(PT_max_best <- results$PT_max[which(results$PT_rmse==min(results$PT_rmse))])
(Penman_max_best <- results$Penman_max[which(results$Penman_rmse==min(results$Penman_rmse))])
(PMKc_max_best <- results$PMKc_max[which(results$PMKc_rmse==min(results$PMKc_rmse))])
(PMgs_max_best <- results$PMgs_max[which(results$PMgs_rmse==min(results$PMgs_rmse))])
(PMgs_min_best <- results$PMgs_min[which(results$PMgs_rmse==min(results$PMgs_rmse))])

head(df)

# PT
df$PT_stress_cor <- ifelse(df$LE_CORR/df$PET_PT < 0, 0, ifelse(df$LE_CORR/df$PET_PT > PT_max_best, PT_max_best, df$LE_CORR/df$PET_PT))

# Penman
df$Penman_stress_cor <- ifelse(df$LE_CORR/df$PET_Penman < 0, 0, ifelse(df$LE_CORR/df$PET_Penman > Penman_max_best, Penman_max_best, df$LE_CORR/df$PET_Penman))

# PM-Kc
df$PM_stress_cor <- ifelse(df$LE_CORR/df$PET_PM < 0, 0, ifelse(df$LE_CORR/df$PET_PM > PMKc_max_best, PMKc_max_best, df$LE_CORR/df$PET_PM))

#PM-gs
df <- df %>%
	mutate(Gs_cor = ifelse(!is.na(Gs_cor_raw), Gs_cor_raw,
												 ifelse(LE_CORR <= 0, PMgs_min_best,
												 			 ifelse(is.infinite(Gs_cor_raw),PMgs_min_best,PMgs_max_best))),
				 Gs_log_cor = log(Gs_cor)
	)

head(df)
## save data
write.csv(df, "FLX_RS_merged_filtered2.csv",row.names = F)

