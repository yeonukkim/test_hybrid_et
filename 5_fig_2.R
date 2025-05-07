rm(list = ls())

# load libraries
library(tidyverse)
library(caret)
library(bigleaf)
library(ggpubr)
library(cowplot)
library(ggpointdensity)
library(viridis)
library(Metrics)
library(ggpmisc)
library(corrplot)
library(hydroGOF)
library(RColorBrewer)
library(sensitivity)
library(boot)
library(ggbreak)

# load dataset
df <- read.csv("FLX_RS_merged_filtered2.csv") 
cp <- bigleaf.constants()$cp

###############################
############ train/ test set
# train/test set index
index <- readRDS("data_partition_index.rds") # trainset
train_set <- df[index,]
test_set <- df[-index,]


########################
##### hybrid model functions
########################
### models
# Model calculation functions
calculate_LE_CORR_hybrid_PT <- function(ML_parameter) {
	train_set %>%
		mutate(LE_CORR_hybrid_PT = PET_PT * ML_parameter) %>%
		pull(LE_CORR_hybrid_PT)
}
calculate_LE_CORR_hybrid_Penman <- function(ML_parameter) {
	train_set %>%
		mutate(LE_CORR_hybrid_Penman = PET_Penman * ML_parameter) %>%
		pull(LE_CORR_hybrid_Penman)
}
calculate_LE_CORR_hybrid_PMkc <- function(ML_parameter) {
	train_set %>%
		mutate(LE_CORR_hybrid_PMkc = PET_PM * ML_parameter) %>%
		pull(LE_CORR_hybrid_PMkc)
}
calculate_LE_CORR_hybrid_PMgs <- function(ML_parameter) {
	train_set %>%
		mutate(LE_CORR_hybrid_PMgs = (Delta * (NETRAD - G_F_MDS) + density * cp * VPD_F / 10 / ra_H) /
					 	(Delta + gamma + gamma /exp(ML_parameter) / ra_H)) %>%
		pull(LE_CORR_hybrid_PMgs)
}
calculate_LE_CORR_hybrid_PMRH_Frh <- function(ML_parameter) {
	train_set %>%
		mutate(LE_CORR_hybrid_PMRH_Frh = rha * Delta * (NETRAD - G_F_MDS) / (rha * Delta + gamma) +
					 	gamma / (rha * Delta + gamma) * ML_parameter * density * cp * esat_air / gamma) %>%
		pull(LE_CORR_hybrid_PMRH_Frh)
}
calculate_LE_CORR_hybrid_PMRH_drh <- function(ML_parameter) {
	train_set %>%
		mutate(LE_CORR_hybrid_PMRH_drh = rha * Delta * (NETRAD - G_F_MDS) / (rha * Delta + gamma) +
					 	gamma / (rha * Delta + gamma) * ML_parameter / ra_H * density * cp * esat_air / gamma) %>%
		pull(LE_CORR_hybrid_PMRH_drh)
}


####################################################
#### training set results with different ntree
####################################################

ntree_values <- c(100,500)

# Initialize the results data frame
results <- data.frame(ntree = ntree_values, Pure_ML = NA, H_PT=NA, H_Penman=NA,
											H_PM_Kc = NA, H_PM_gs = NA, H_PMRH_Frh = NA ,H_PMRH_drh = NA)

# Loop over the ntree values
for (i in seq_along(ntree_values)) {
  ntree <- ntree_values[i]
	# Read the model from the .rds file
	RF_LE_CORR <- readRDS(paste0("./ntree_models/RF_LE_CORR_ntree_", ntree, ".rds"))
	RF_PT_stress_cor <- readRDS(paste0("./ntree_models/RF_PT_stress_cor_ntree_", ntree, ".rds"))
	RF_Penman_stress_cor <- readRDS(paste0("./ntree_models/RF_Penman_stress_cor_ntree_", ntree, ".rds"))
	RF_PM_stress_cor <- readRDS(paste0("./ntree_models/RF_PM_stress_cor_ntree_", ntree, ".rds"))
	RF_Gs_log_cor <- readRDS(paste0("./ntree_models/RF_Gs_log_cor_ntree_", ntree, ".rds"))
	RF_Frh_cor <- readRDS(paste0("./ntree_models/RF_Frh_cor_ntree_", ntree, ".rds"))
	RF_drh_cor <- readRDS(paste0("./ntree_models/RF_drh_cor_ntree_", ntree, ".rds"))
	# Make predictions on the train set
	LE_CORR_RF <- predict(RF_LE_CORR, train_set, na.action = na.pass)
	PT_stress_cor_RF <- predict(RF_PT_stress_cor, train_set, na.action = na.pass)
	Penman_stress_cor_RF <- predict(RF_Penman_stress_cor, train_set, na.action = na.pass)
	PM_stress_cor_RF <- predict(RF_PM_stress_cor, train_set, na.action = na.pass)
	Gs_log_cor_RF <- predict(RF_Gs_log_cor, train_set, na.action = na.pass)
	Frh_cor_RF <- predict(RF_Frh_cor, train_set, na.action = na.pass)
	drh_cor_RF <- predict(RF_drh_cor, train_set, na.action = na.pass)
	# hybrid model application
	LE_CORR_hybrid_PT <- calculate_LE_CORR_hybrid_PT(PT_stress_cor_RF)
	LE_CORR_hybrid_Penman <- calculate_LE_CORR_hybrid_Penman(Penman_stress_cor_RF)
	LE_CORR_hybrid_PMkc <- calculate_LE_CORR_hybrid_PMkc(PM_stress_cor_RF)
	LE_CORR_hybrid_PMgs <- calculate_LE_CORR_hybrid_PMgs(Gs_log_cor_RF)
	LE_CORR_hybrid_PMRH_Frh <- calculate_LE_CORR_hybrid_PMRH_Frh(Frh_cor_RF)
	LE_CORR_hybrid_PMRH_drh <- calculate_LE_CORR_hybrid_PMRH_drh(drh_cor_RF)
	# Calculate RMSE and store it in the results data frame
	results$Pure_ML[i] <- rmse(LE_CORR_RF, train_set$LE_CORR)
	results$H_PT[i] <- rmse(LE_CORR_hybrid_PT, train_set$LE_CORR)
	results$H_Penman[i] <- rmse(LE_CORR_hybrid_Penman, train_set$LE_CORR)
	results$H_PM_Kc[i] <- rmse(LE_CORR_hybrid_PMkc, train_set$LE_CORR)
	results$H_PM_gs[i] <- rmse(LE_CORR_hybrid_PMgs, train_set$LE_CORR)
	results$H_PMRH_Frh[i] <- rmse(LE_CORR_hybrid_PMRH_Frh, train_set$LE_CORR)
	results$H_PMRH_drh[i] <- rmse(LE_CORR_hybrid_PMRH_drh, train_set$LE_CORR)
}


####################################################
#### training set results with different mtry 
####################################################

# Read the pure ML model (including mtry tuning results)
RF_LE_CORR <- readRDS("./ntree_models/RF_LE_CORR_ntree_100.rds")

improvements <- data.frame(
	ImprovementType =c('a','b','c'),
	ImprovementValue = c(
		# model improvement due to increase ntree 100 -> 500
		results$Pure_ML[which(results$ntree == 500)] - results$Pure_ML[which(results$ntree == 100)],
		# model improvement due to mtry tuning (default: 9/3 = 3):
		min(RF_LE_CORR$results$RMSE) - RF_LE_CORR$results$RMSE[3], 
		# model improvement using hybrid model
		results$H_PMRH_Frh[1] - results$Pure_ML[1]
	)
)

p2<- ggplot(improvements, aes(x = ImprovementType, y = ImprovementValue, fill = ImprovementType)) +
	geom_bar(stat = "identity",width = 0.5) +
	labs(
		y = expression("Decrease in RMSE (W " * m^-2 * ")")
	) +
	scale_fill_manual(values = c("#009E73", "#F0E442", "black")) +
	theme_bw() +
	theme(
		axis.text.x = element_text(angle = 0, size = 12, vjust = 0.5),
		legend.position = "none",
		axis.title.x.bottom = element_blank()
	) +
	geom_hline(yintercept = 0) +
	geom_vline(xintercept = 2.5) +
	scale_x_discrete(labels = c(
		expression(atop("Increase ntree", "(100" %->% "500)")),
		expression(atop("Tuning mtry", "(default" %->% "best)")),
		expression(atop("Employing hybrid", "(Pure-ML" %->% "H-" * PM[RH] * "-" * F[RH] * ")"))
	))


tiff(filename = "./Figures/tif_Figure2.tif",width = 7, height = 4,units = "in", res = 300)
p2
dev.off()


