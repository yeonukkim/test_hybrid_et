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

# load dataset
df <- read.csv("FLX_RS_merged_filtered2.csv") 
cp <- bigleaf.constants()$cp

#######################
#### model results
#######################
RF_LE_CORR <- readRDS("./ntree_models/RF_LE_CORR_ntree_500.rds")
df$LE_CORR_RF <- predict(RF_LE_CORR, df, na.action = na.pass)

RF_PT_stress_cor <- readRDS("./ntree_models/RF_PT_stress_cor_ntree_500.rds")
df$PT_stress_cor_RF <- predict(RF_PT_stress_cor, df, na.action = na.pass)

RF_Penman_stress_cor <- readRDS("./ntree_models/RF_Penman_stress_cor_ntree_500.rds")
df$Penman_stress_cor_RF <- predict(RF_Penman_stress_cor, df, na.action = na.pass)

RF_PM_stress_cor <- readRDS("./ntree_models/RF_PM_stress_cor_ntree_500.rds")
df$PM_stress_cor_RF <- predict(RF_PM_stress_cor, df, na.action = na.pass)

RF_Gs_log_cor <- readRDS("./ntree_models/RF_Gs_log_cor_ntree_500.rds")
df$Gs_log_cor_RF <- predict(RF_Gs_log_cor, df, na.action = na.pass)

RF_Frh_cor <- readRDS("./ntree_models/RF_Frh_cor_ntree_500.rds")
df$Frh_cor_RF <- predict(RF_Frh_cor, df, na.action = na.pass)

RF_drh_cor <- readRDS("./ntree_models/RF_drh_cor_ntree_500.rds")
df$drh_cor_RF <- predict(RF_drh_cor, df, na.action = na.pass)


########################
## hybrid LE
########################

### models
# Model calculation functions
calculate_LE_CORR_hybrid_PT <- function(ML_parameter) {
	df %>%
		mutate(LE_CORR_hybrid_PT = PET_PT * ML_parameter) %>%
		pull(LE_CORR_hybrid_PT)
}

calculate_LE_CORR_hybrid_Penman <- function(ML_parameter) {
	df %>%
		mutate(LE_CORR_hybrid_Penman = PET_Penman * ML_parameter) %>%
		pull(LE_CORR_hybrid_Penman)
}

calculate_LE_CORR_hybrid_PMkc <- function(ML_parameter) {
	df %>%
		mutate(LE_CORR_hybrid_PMkc = PET_PM * ML_parameter) %>%
		pull(LE_CORR_hybrid_PMkc)
}

calculate_LE_CORR_hybrid_PMgs <- function(ML_parameter) {
	df %>%
		mutate(LE_CORR_hybrid_PMgs = (Delta * (NETRAD - G_F_MDS) + density * cp * VPD_F / 10 / ra_H) / 
					 	(Delta + gamma + gamma /exp(ML_parameter) / ra_H)) %>%
		pull(LE_CORR_hybrid_PMgs)
}

calculate_LE_CORR_hybrid_PMRH_Frh <- function(ML_parameter) {
	df %>%
		mutate(LE_CORR_hybrid_PMRH_Frh = rha * Delta * (NETRAD - G_F_MDS) / (rha * Delta + gamma) +
					 	gamma / (rha * Delta + gamma) * ML_parameter * density * cp * esat_air / gamma) %>%
		pull(LE_CORR_hybrid_PMRH_Frh)
}

calculate_LE_CORR_hybrid_PMRH_drh <- function(ML_parameter) {
	df %>%
		mutate(LE_CORR_hybrid_PMRH_drh = rha * Delta * (NETRAD - G_F_MDS) / (rha * Delta + gamma) +
					 	gamma / (rha * Delta + gamma) * ML_parameter / ra_H * density * cp * esat_air / gamma) %>%
		pull(LE_CORR_hybrid_PMRH_drh)
}


# apply
df$LE_CORR_hybrid_PT <- calculate_LE_CORR_hybrid_PT(df$PT_stress_cor_RF)
df$LE_CORR_hybrid_Penman <- calculate_LE_CORR_hybrid_Penman(df$Penman_stress_cor_RF)
df$LE_CORR_hybrid_PMkc <- calculate_LE_CORR_hybrid_PMkc(df$PM_stress_cor_RF)
df$LE_CORR_hybrid_PMgs <- calculate_LE_CORR_hybrid_PMgs(df$Gs_log_cor_RF)
df$LE_CORR_hybrid_PMRH_Frh <- calculate_LE_CORR_hybrid_PMRH_Frh(df$Frh_cor_RF)
df$LE_CORR_hybrid_PMRH_drh <- calculate_LE_CORR_hybrid_PMRH_drh(df$drh_cor_RF)

###############################
############ train/ test set
# train/test set index
index <- readRDS("data_partition_index.rds") # trainset
train_set <- df[index,]
test_set <- df[-index,]

######################################
####### test set performance
######################################

a <- test_set %>% 
	ggplot(aes(LE_CORR_RF, LE_CORR)) +
	geom_point(alpha = 0.1) +
	annotate("text", x = -Inf, y = 320, hjust = -0.1, size = 3,
					 label = sprintf("RMSE = %.2f", rmse(test_set$LE_CORR, test_set$LE_CORR_RF))) +
	annotate("text", x = -Inf, y = 280, hjust = -0.1, size = 3,
					 label = sprintf("sigma[m]/sigma[o] == %.3f", sd(test_set$LE_CORR_RF)/sd(test_set$LE_CORR)), parse = TRUE) +
	geom_abline(lty = 2) +
	# geom_smooth(method = "lm", color = "orange") +
	ylab(expression(LE[obs]*" (W "*m^-2*")")) +
	xlab(expression(LE[model]*" (W "*m^-2*")")) +
	ggtitle(expression("a. "*italic("Pure-ML"))) +
	theme_bw() + xlim(c(-50,300)) + ylim(c(-50,350)) +
	theme(title = element_text(size = 8),
				axis.text = element_text(size = 7),
				axis.title = element_text(size = 7))


b <- test_set %>% 
	ggplot(aes(LE_CORR_hybrid_Penman, LE_CORR)) +
	geom_point(alpha = 0.1) +
	annotate("text", x = -Inf, y = 320, hjust = -0.1, size = 3,
					 label = sprintf("RMSE = %.2f", rmse(test_set$LE_CORR, test_set$LE_CORR_hybrid_Penman))) +
	annotate("text", x = -Inf, y = 280, hjust = -0.1, size = 3,
					 label = sprintf("sigma[m]/sigma[o] == %.3f", sd(test_set$LE_CORR_hybrid_Penman)/sd(test_set$LE_CORR)), 
					 parse = TRUE) +
	geom_abline(lty = 2) +
	# geom_smooth(method = "lm", color = "orange") +
	ylab(expression(LE[obs]*" (W "*m^-2*")")) +
	xlab(expression(LE[model]*" (W "*m^-2*")")) +
	ggtitle(expression("b. "*italic("H-Penman"))) +
	theme_bw() + xlim(c(-50,300)) + ylim(c(-50,350))+
	theme(title = element_text(size = 8),
				axis.text = element_text(size = 7),
				axis.title = element_text(size = 7))


c <- test_set %>% 
	ggplot(aes(LE_CORR_hybrid_PT, LE_CORR)) +
	geom_point(alpha = 0.1) +
	annotate("text", x = -Inf, y = 320, hjust = -0.1, size = 3,
					 label = sprintf("RMSE = %.2f", rmse(test_set$LE_CORR, test_set$LE_CORR_hybrid_PT))) +
	annotate("text", x = -Inf, y = 280, hjust = -0.1, size = 3,
					 label = sprintf("sigma[m]/sigma[o] == %.3f", sd(test_set$LE_CORR_hybrid_PT)/sd(test_set$LE_CORR)), parse = TRUE) +
	geom_abline(lty = 2) +
	# geom_smooth(method = "lm", color = "orange") +
	ylab(expression(LE[obs]*" (W "*m^-2*")")) +
	xlab(expression(LE[model]*" (W "*m^-2*")")) +
	ggtitle(expression("c. "*italic("H-PT"))) +
	theme_bw() + xlim(c(-50,300)) + ylim(c(-50,350))+
	theme(title = element_text(size = 8),
				axis.text = element_text(size = 7),
				axis.title = element_text(size = 7))


d <- test_set %>% 
	ggplot(aes(LE_CORR_hybrid_PMgs, LE_CORR)) +
	geom_point(alpha = 0.1) +
	annotate("text", x = -Inf, y = 320, hjust = -0.1, size = 3,
					 label = sprintf("RMSE = %.2f", rmse(test_set$LE_CORR, test_set$LE_CORR_hybrid_PMgs))) +
	annotate("text", x = -Inf, y = 280, hjust = -0.1, size = 3,
					 label = sprintf("sigma[m]/sigma[o] == %.3f", sd(test_set$LE_CORR_hybrid_PMgs)/sd(test_set$LE_CORR)), parse = TRUE) +
	geom_abline(lty = 2) +
	# geom_smooth(method = "lm", color = "orange") +
	ylab(expression(LE[obs]*" (W "*m^-2*")")) +
	xlab(expression(LE[model]*" (W "*m^-2*")")) +
	ggtitle(expression("d. "*italic("H-PM-"*g[s]))) +
	theme_bw() + xlim(c(-50,300)) + ylim(c(-50,350))+
	theme(title = element_text(size = 8),
				axis.text = element_text(size = 7),
				axis.title = element_text(size = 7))

e <- test_set %>% 
	ggplot(aes(LE_CORR_hybrid_PMkc, LE_CORR)) +
	geom_point(alpha = 0.1) +
	annotate("text", x = -Inf, y = 320, hjust = -0.1, size = 3,
					 label = sprintf("RMSE = %.2f", rmse(test_set$LE_CORR, test_set$LE_CORR_hybrid_PMkc))) +
	annotate("text", x = -Inf, y = 280, hjust = -0.1, size = 3,
					 label = sprintf("sigma[m]/sigma[o] == %.3f", sd(test_set$LE_CORR_hybrid_PMkc)/sd(test_set$LE_CORR)), parse = TRUE) +
	geom_abline(lty = 2) +
	# geom_smooth(method = "lm", color = "orange") +
	ylab(expression(LE[obs]*" (W "*m^-2*")")) +
	xlab(expression(LE[model]*" (W "*m^-2*")")) +
	ggtitle(expression("e. "*italic("H-PM-"*K[c]))) +
	theme_bw() + xlim(c(-50,300)) + ylim(c(-50,350))+
	theme(title = element_text(size = 8),
				axis.text = element_text(size = 7),
				axis.title = element_text(size = 7))


f <- test_set %>% 
	ggplot(aes(LE_CORR_hybrid_PMRH_drh, LE_CORR)) +
	geom_point(alpha = 0.1) +
	annotate("text", x = -Inf, y = 320, hjust = -0.1, size = 3,
					 label = sprintf("RMSE = %.2f", rmse(test_set$LE_CORR, test_set$LE_CORR_hybrid_PMRH_drh))) +
	annotate("text", x = -Inf, y = 280, hjust = -0.1, size = 3,
					 label = sprintf("         = %.3f", sd(test_set$LE_CORR_hybrid_PMRH_drh)/sd(test_set$LE_CORR))) +
	annotate("text", x = -Inf, y = 280, hjust = -0.1, size = 3,
					 label = sprintf("sigma[m]/sigma[o]"), parse = TRUE) +
	geom_abline(lty = 2) +
	# geom_smooth(method = "lm", color = "orange") +
	ylab(expression(LE[obs]*" (W "*m^-2*")")) +
	xlab(expression(LE[model]*" (W "*m^-2*")")) +
	ggtitle(expression("f. "*italic("H-"*PM[RH]*'-'*d[RH]))) +
	theme_bw() + xlim(c(-50,300)) + ylim(c(-50,350))+
	theme(title = element_text(size = 8),
				axis.text = element_text(size = 7),
				axis.title = element_text(size = 7))


g <- test_set %>% 
	ggplot(aes(LE_CORR_hybrid_PMRH_Frh, LE_CORR)) +
	geom_point(alpha = 0.1) +
	annotate("text", x = -Inf, y = 320, hjust = -0.1, size = 3,
					 label = sprintf("RMSE = %.2f", rmse(test_set$LE_CORR, test_set$LE_CORR_hybrid_PMRH_Frh))) +
	annotate("text", x = -Inf, y = 280, hjust = -0.1, size = 3,
					 label = sprintf("         = %.3f", sd(test_set$LE_CORR_hybrid_PMRH_Frh)/sd(test_set$LE_CORR))) +
	annotate("text", x = -Inf, y = 280, hjust = -0.1, size = 3,
					 label = sprintf("sigma[m]/sigma[o]"), parse = TRUE) +
	geom_abline(lty = 2) +
	# geom_smooth(method = "lm", color = "orange") +
	ylab(expression(LE[obs]*" (W "*m^-2*")")) +
	xlab(expression(LE[model]*" (W "*m^-2*")")) +
	ggtitle(expression("g. "*italic("H-"*PM[RH]*'-'*F[RH]))) +
	theme_bw() + xlim(c(-50,300)) + ylim(c(-50,350))+
	theme(title = element_text(size = 8),
				axis.text = element_text(size = 7),
				axis.title = element_text(size = 7))


figure_a <- plot_grid(plot_grid(NULL,a,b,c,NULL,ncol = 5, rel_widths = c(0.5,1,1,1,0.5)),
											plot_grid(d,e,f,g,ncol = 4),ncol = 1)
	
tiff(filename = "./Figures/tif_Figure3.tif",width = 7.5, height = 4,units = "in", res = 300)
figure_a
dev.off()


###############################################
### performance in different IGBP (figure 5)
###############################################

meta <- read.csv('IGBP_from_AA.csv')
meta$IGBP2 <- ifelse(meta$IGBP %in% c("OSH","CSH"),"OSH/CSH (n = 7)",
										 ifelse(meta$IGBP == "CRO", "CRO (n = 7)",
										 			 ifelse(meta$IGBP %in% c("DBF","ENF"), "DBF/ENF (n = 11)",
										 			 			 			 ifelse(meta$IGBP == "GRA", "GRA (n = 8)", "WET (n = 7)"
										 			 			 			 			 ))))

test_set2 <- left_join(test_set, meta, by='site')

aa <-  test_set2 %>% 
	summarise(H1 = rmse(LE_CORR,LE_CORR_hybrid_Penman)- rmse(LE_CORR,LE_CORR_RF),
						H2 = rmse(LE_CORR,LE_CORR_hybrid_PT)- rmse(LE_CORR,LE_CORR_RF),
						H3 = rmse(LE_CORR,LE_CORR_hybrid_PMgs)- rmse(LE_CORR,LE_CORR_RF),
						H4 = rmse(LE_CORR,LE_CORR_hybrid_PMkc)- rmse(LE_CORR,LE_CORR_RF),
						H5 = rmse(LE_CORR,LE_CORR_hybrid_PMRH_drh)- rmse(LE_CORR,LE_CORR_RF),
						H6 = rmse(LE_CORR,LE_CORR_hybrid_PMRH_Frh)- rmse(LE_CORR,LE_CORR_RF)
	)
aa <- aa %>% gather(key = "model",value = "rmse")
aa$IGBP2 <- "All (n = 40)"


bb <-  test_set2 %>% 
	group_by(IGBP2) %>%
	summarise(H1 = rmse(LE_CORR,LE_CORR_hybrid_Penman)- rmse(LE_CORR,LE_CORR_RF),
						H2 = rmse(LE_CORR,LE_CORR_hybrid_PT)- rmse(LE_CORR,LE_CORR_RF),
						H3 = rmse(LE_CORR,LE_CORR_hybrid_PMgs)- rmse(LE_CORR,LE_CORR_RF),
						H4 = rmse(LE_CORR,LE_CORR_hybrid_PMkc)- rmse(LE_CORR,LE_CORR_RF),
						H5 = rmse(LE_CORR,LE_CORR_hybrid_PMRH_drh)- rmse(LE_CORR,LE_CORR_RF),
						H6 = rmse(LE_CORR,LE_CORR_hybrid_PMRH_Frh)- rmse(LE_CORR,LE_CORR_RF)
	)
bb <- bb %>% gather(key = "model",value = "rmse", - c("IGBP2"))


use <- test_set2 %>% 
	group_by(IGBP2,site) %>%
	summarise(H1 = rmse(LE_CORR,LE_CORR_hybrid_Penman)- rmse(LE_CORR,LE_CORR_RF),
						H2 = rmse(LE_CORR,LE_CORR_hybrid_PT)- rmse(LE_CORR,LE_CORR_RF),
						H3 = rmse(LE_CORR,LE_CORR_hybrid_PMgs)- rmse(LE_CORR,LE_CORR_RF),
						H4 = rmse(LE_CORR,LE_CORR_hybrid_PMkc)- rmse(LE_CORR,LE_CORR_RF),
						H5 = rmse(LE_CORR,LE_CORR_hybrid_PMRH_drh)- rmse(LE_CORR,LE_CORR_RF),
						H6 = rmse(LE_CORR,LE_CORR_hybrid_PMRH_Frh)- rmse(LE_CORR,LE_CORR_RF)
						) %>%
	ungroup()

use_long <- use %>%
	gather(key = "model",value = "rmse", -c("IGBP2","site")) 

use_long2 <- use_long %>%
	mutate(IGBP2 = "All (n = 40)")

use_long <- rbind(use_long,use_long2)

p <- use_long %>%
	ggplot(aes(model,rmse,fill = model)) +
	geom_hline(yintercept = 0) +
	geom_boxplot() +
	facet_wrap(vars(IGBP2), scales = "free_y") +
	theme_bw() +
	theme(axis.title.x = element_blank(),
				legend.text = element_text(hjust = 0, size = 8),
				legend.position = "bottom") +
	ylab(expression("RMSE difference (Hybrid - pure ML) "*" (W "*m^-2*")")) +
	scale_fill_manual(name = "",
		labels = c("H1: H-Penman", "H2: H-PT", expression("H3: H-PM-"*g[s]),
							 expression("H4: H-PM-"*K[c]), expression("H5: H-"*PM[RH]*"-"*d[RH]), 
							 expression("H6: H-"*PM[RH]*"-"*F[RH])
		),
		values = c("#999999","#E69F00", "#009E73","#0072B2", "#CC79A7", "#D55E00")
	) +
	geom_point(data = aa, fill = "white", size = 2 , shape = 24, color = "black") +
	geom_point(data = bb, fill = "white", size = 2 , shape = 24, color = "black") +
	guides(fill = guide_legend(nrow = 1))


tiff(filename = "./Figures/tif_Figure4.tif",width = 7.5, height = 4.5,units = "in", res = 300)
p
dev.off()

######################################
### sensitivity analysis
######################################

anal <- data.frame(models = c("Pure-ML","H-PT","H-Penman","H-PM-Kc","H-PM-gs","H-PMRH-Frh","H-PMRH-drh"),
									 partial_stdev = c(1,1,1,1,1,1,1),
									 stdev_ratio = c(sd(test_set$LE_CORR_RF)/sd(test_set$LE_CORR),
									 					sd(test_set$LE_CORR_hybrid_PT)/sd(test_set$LE_CORR),
									 					sd(test_set$LE_CORR_hybrid_Penman)/sd(test_set$LE_CORR),
									 					sd(test_set$LE_CORR_hybrid_PMkc)/sd(test_set$LE_CORR),
									 					sd(test_set$LE_CORR_hybrid_PMgs)/sd(test_set$LE_CORR),
									 					sd(test_set$LE_CORR_hybrid_PMRH_Frh)/sd(test_set$LE_CORR),
									 					sd(test_set$LE_CORR_hybrid_PMRH_drh)/sd(test_set$LE_CORR)
									 ),
									 mean_bias = c(bias(test_set$LE_CORR_RF, test_set$LE_CORR),
									 							bias(test_set$LE_CORR_hybrid_PT, test_set$LE_CORR),
									 							bias(test_set$LE_CORR_hybrid_Penman, test_set$LE_CORR),
									 							bias(test_set$LE_CORR_hybrid_PMkc, test_set$LE_CORR),
									 							bias(test_set$LE_CORR_hybrid_PMgs, test_set$LE_CORR),
									 							bias(test_set$LE_CORR_hybrid_PMRH_Frh, test_set$LE_CORR),
									 							bias(test_set$LE_CORR_hybrid_PMRH_drh, test_set$LE_CORR)
									 ),
									 rmse = c(rmse(test_set$LE_CORR,test_set$LE_CORR_RF),
									 				 rmse(test_set$LE_CORR,test_set$LE_CORR_hybrid_PT),
									 				 rmse(test_set$LE_CORR,test_set$LE_CORR_hybrid_Penman),
									 				 rmse(test_set$LE_CORR,test_set$LE_CORR_hybrid_PMkc),
									 				 rmse(test_set$LE_CORR,test_set$LE_CORR_hybrid_PMgs),
									 				 rmse(test_set$LE_CORR,test_set$LE_CORR_hybrid_PMRH_Frh),
									 				 rmse(test_set$LE_CORR,test_set$LE_CORR_hybrid_PMRH_drh)
									 ),
									 parameter_R2 = c(cor(test_set$LE_CORR,test_set$LE_CORR_RF)^2,
									 								 cor(test_set$PT_stress_cor,test_set$PT_stress_cor_RF)^2,
									 								 cor(test_set$Penman_stress_cor,test_set$Penman_stress_cor_RF)^2,
									 								 cor(test_set$PM_stress_cor,test_set$PM_stress_cor_RF)^2,
									 								 cor(test_set$Gs_log_cor,test_set$Gs_log_cor_RF)^2,
									 								 cor(test_set$Frh_cor,test_set$Frh_cor_RF)^2,
									 								 cor(test_set$drh_cor,test_set$drh_cor_RF)^2
									 )
)


anal$partial_stdev[1] <- sd(test_set$LE_CORR_RF)

#PT
anal$partial_stdev[2] <- sd(test_set$PT_stress_cor_RF) *
	test_set %>%
	summarise(sens = sqrt(mean((PET_PT)^2))
						) %>%
	pull(sens)

#Penman
anal$partial_stdev[3] <- sd(test_set$Penman_stress_cor_RF) *
	test_set %>%
	summarise(sens = sqrt(mean((PET_Penman)^2))) %>%
	pull(sens)

#PMKc
anal$partial_stdev[4] <- sd(test_set$PM_stress_cor_RF)*
	test_set %>%
	summarise(sens = sqrt(mean((PET_PM)^2))) %>%
	pull(sens)

# gs
anal$partial_stdev[5] <- sd(test_set$Gs_log_cor_RF)*
	test_set %>%
	summarise(sens = sqrt(mean(((Delta*AE+density*cp*VPD_F/10/ra_H)*gamma/ra_H/
													(Delta+gamma+gamma/exp(Gs_log_cor_RF)/ra_H)^2/(exp(Gs_log_cor_RF)))^2
													)
												)
	) %>%
	pull(sens)
	

#Frh
anal$partial_stdev[6] <- 	sd(test_set$Frh_cor_RF)*
	test_set %>%
	summarise(sens = sqrt(mean((gamma/(rha * Delta + gamma) * density * cp * esat_air / gamma))^2
												)
						) %>%
	pull(sens)

#drh
anal$partial_stdev[7] <- sd(test_set$drh_cor_RF)*
	test_set %>%
	summarise(sens = sqrt(mean((gamma/(rha * Delta + gamma) * density * cp * esat_air / gamma/ra_H))^2
												)
						) %>%
	pull(sens)


############# results


p0 <- anal %>%
	filter(models != "Pure-ML") %>%
	ggplot(aes(parameter_R2, rmse)) +
	geom_point(size = 1) +
	geom_text(aes(label = case_when(
		models == "Pure-ML" ~ "italic('Pure-ML ')",
		models == "H-Penman" ~ "italic('H-Penman ')",
		models == "H-PT" ~ "italic('H-PT  ')",
		models == "H-PM-Kc" ~ "italic('H-PM-' * K[c])",
		models == "H-PM-gs" ~ "italic('H-PM-' * g[s])",
		models == "H-PMRH-drh" ~ "italic('H-'* PM[RH]*'-'*d[RH])",
		models == "H-PMRH-Frh" ~ "italic('H-'* PM[RH]*'-'*F[RH])"
	)), parse = TRUE, hjust = 0.5, vjust = 1.5, angle = 0,size = 3) +
	theme_bw() +
	ylab(expression("RMSE "*" (W " * m^-2 * ")")) +
	xlab(expression("ML estimated intermediate parameter "*R^2)) +
	stat_cor() +
	xlim(c(0.42,0.72)) + ylim(c(17.6,19.25))

tiff(filename = "./Figures/tif_Figure5.tif",width = 4, height = 3.5,units = "in", res = 300)
p0
dev.off()


p <- ggplot(anal, aes(partial_stdev, rmse)) +
	geom_point(size = 1) +
	geom_text(aes(label = case_when(
		models == "Pure-ML" ~ "italic('Pure-ML ')",
		models == "H-Penman" ~ "italic('H-Penman ')",
		models == "H-PT" ~ "italic('H-PT  ')",
		models == "H-PM-Kc" ~ "italic('H-PM-' * K[c])",
		models == "H-PM-gs" ~ "italic('H-PM-' * g[s])",
		models == "H-PMRH-drh" ~ "italic('H-'* PM[RH]*'-'*d[RH])",
		models == "H-PMRH-Frh" ~ "italic('H-'* PM[RH]*'-'*F[RH])"
	)), parse = TRUE, hjust = 0.5, vjust = 1.5, angle = 20, size = 3) +
	theme_bw() +
	ylab(expression("RMSE "*" (W " * m^-2 * ")")) +
	xlab(expression("Global Sensitivity Index" * " (W " * m^-2 * ")")) +
	geom_hline(yintercept = anal$rmse[1], lty = 3) +
	geom_vline(xintercept = anal$partial_stdev[1], lty = 3) +
	stat_cor() +
	xlim(c(12,47)) + ylim(c(17.60,19.3)) +
	ggtitle(expression("GSI vs. ET models' RMSE"))

p2 <- ggplot(anal, aes(partial_stdev, stdev_ratio)) +
	geom_point(size = 1) +
	geom_text(aes(label = case_when(
		models == "Pure-ML" ~ "italic('Pure-ML ')",
		models == "H-Penman" ~ "italic('H-Penman ')",
		models == "H-PT" ~ "italic('H-PT  ')",
		models == "H-PM-Kc" ~ "italic('H-PM-' * K[c])",
		models == "H-PM-gs" ~ "italic('H-PM-' * g[s])",
		models == "H-PMRH-drh" ~ "italic('H-'* PM[RH]*'-'*d[RH])",
		models == "H-PMRH-Frh" ~ "italic('H-'* PM[RH]*'-'*F[RH]*' ')"
	)), parse = TRUE, hjust = 0.5, vjust = 1.5, angle = 20, size = 3) +
	theme_bw() +
	ylab(expression(italic(sigma[m]/sigma[o]))) +
	xlab(expression("Global Sensitivity Index" * " (W " * m^-2 * ")")) +
	geom_hline(yintercept = anal$stdev_ratio[1], lty = 3) +
	geom_vline(xintercept = anal$partial_stdev[1], lty = 3) +
	stat_cor() +
	xlim(c(12,47)) + ylim(c(0.87,0.95))+
	ggtitle(expression("GSI vs. ET models' "*italic(sigma[m]/sigma[o])))



figure_bc <- plot_grid(p,p2,ncol = 2,labels = "auto")


tiff(filename = "./Figures/tif_Figure6.tif",width = 7, height = 3,units = "in", res = 300)
figure_bc
dev.off()




#################################
#### extreme condition analysis: all variables
################################
#### 1%
extreme <- df %>% 
	group_by(site)%>%
	summarise(VWC_03 = quantile(soil_moisture,0.01),
						RH_03 = quantile(rha,0.01),
						T_03 = quantile(TA_F,0.01),
						FPAR_03 = quantile(Fpar,0.01),
						SW_IN_03 = quantile(SW_IN_F_MDS,0.01),
						ga_H_03 = quantile(ga_H,0.01),
						VegWC_03 = quantile(VegWC, 0.01),
						CO2_03 = quantile(CO2_F_MDS, 0.01),
						
						VWC_97 = quantile(soil_moisture,0.99),
						RH_97 = quantile(rha,0.99),
						T_97 = quantile(TA_F,0.99),
						FPAR_97 = quantile(Fpar,0.99),
						SW_IN_97 = quantile(SW_IN_F_MDS,0.99),
						ga_H_97 = quantile(ga_H,0.99),
						VegWC_97 = quantile(VegWC, 0.99),
						CO2_97 = quantile(CO2_F_MDS, 0.99)
	)	

df2 <- df %>% 
	left_join(extreme)

test_set2 <- df2[-index,]

# soil moisture
aa <- test_set2 %>% filter(soil_moisture < VWC_03) %>%
	filter(site != "US-CS2") %>%
	summarise(H1 =  rmse( LE_CORR,LE_CORR_hybrid_Penman) -  rmse( LE_CORR,LE_CORR_RF),
						H2 =  rmse( LE_CORR,LE_CORR_hybrid_PT)-  rmse( LE_CORR,LE_CORR_RF),
						H3 =  rmse( LE_CORR,LE_CORR_hybrid_PMgs)-  rmse( LE_CORR,LE_CORR_RF),
						H4 =  rmse( LE_CORR,LE_CORR_hybrid_PMkc)-  rmse( LE_CORR,LE_CORR_RF),
						H5 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_drh)-  rmse( LE_CORR,LE_CORR_RF),
						H6 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_Frh)-  rmse( LE_CORR,LE_CORR_RF)
	)%>%
	gather(key = "model", value = "VWC_03")

aa <- aa %>% left_join(
	test_set2 %>% filter(soil_moisture > VWC_97) %>%
		filter(site != "US-CS2") %>%
		summarise(H1 =  rmse( LE_CORR,LE_CORR_hybrid_Penman) -  rmse( LE_CORR,LE_CORR_RF),
							H2 =  rmse( LE_CORR,LE_CORR_hybrid_PT)-  rmse( LE_CORR,LE_CORR_RF),
							H3 =  rmse( LE_CORR,LE_CORR_hybrid_PMgs)-  rmse( LE_CORR,LE_CORR_RF),
							H4 =  rmse( LE_CORR,LE_CORR_hybrid_PMkc)-  rmse( LE_CORR,LE_CORR_RF),
							H5 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_drh)-  rmse( LE_CORR,LE_CORR_RF),
							H6 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_Frh)-  rmse( LE_CORR,LE_CORR_RF)
		)%>%
		gather(key = "model", value = "VWC_97")
)

#relative humidity
aa <- aa %>% left_join(
	test_set2 %>% filter(rha < RH_03) %>%
		filter(site != "US-CS2") %>%
		summarise(H1 =  rmse( LE_CORR,LE_CORR_hybrid_Penman) -  rmse( LE_CORR,LE_CORR_RF),
							H2 =  rmse( LE_CORR,LE_CORR_hybrid_PT)-  rmse( LE_CORR,LE_CORR_RF),
							H3 =  rmse( LE_CORR,LE_CORR_hybrid_PMgs)-  rmse( LE_CORR,LE_CORR_RF),
							H4 =  rmse( LE_CORR,LE_CORR_hybrid_PMkc)-  rmse( LE_CORR,LE_CORR_RF),
							H5 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_drh)-  rmse( LE_CORR,LE_CORR_RF),
							H6 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_Frh)-  rmse( LE_CORR,LE_CORR_RF)
		)%>%
		gather(key = "model", value = "RH_03")
)

aa <- aa %>% left_join(
	test_set2 %>% filter(rha > RH_97) %>%
		filter(site != "US-CS2") %>%
		summarise(H1 =  rmse( LE_CORR,LE_CORR_hybrid_Penman) -  rmse( LE_CORR,LE_CORR_RF),
							H2 =  rmse( LE_CORR,LE_CORR_hybrid_PT)-  rmse( LE_CORR,LE_CORR_RF),
							H3 =  rmse( LE_CORR,LE_CORR_hybrid_PMgs)-  rmse( LE_CORR,LE_CORR_RF),
							H4 =  rmse( LE_CORR,LE_CORR_hybrid_PMkc)-  rmse( LE_CORR,LE_CORR_RF),
							H5 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_drh)-  rmse( LE_CORR,LE_CORR_RF),
							H6 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_Frh)-  rmse( LE_CORR,LE_CORR_RF)
		)%>%
		gather(key = "model", value = "RH_97")
)

#temperature
aa <- aa %>% left_join(
	test_set2 %>% filter(TA_F < T_03) %>%
		filter(site != "US-CS2") %>%
		summarise(H1 =  rmse( LE_CORR,LE_CORR_hybrid_Penman) -  rmse( LE_CORR,LE_CORR_RF),
							H2 =  rmse( LE_CORR,LE_CORR_hybrid_PT)-  rmse( LE_CORR,LE_CORR_RF),
							H3 =  rmse( LE_CORR,LE_CORR_hybrid_PMgs)-  rmse( LE_CORR,LE_CORR_RF),
							H4 =  rmse( LE_CORR,LE_CORR_hybrid_PMkc)-  rmse( LE_CORR,LE_CORR_RF),
							H5 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_drh)-  rmse( LE_CORR,LE_CORR_RF),
							H6 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_Frh)-  rmse( LE_CORR,LE_CORR_RF)
		)%>%
		gather(key = "model", value = "T_03")
)

aa <- aa %>% left_join(
	test_set2 %>% filter(TA_F > T_97) %>%
		filter(site != "US-CS2") %>%
		summarise(H1 =  rmse( LE_CORR,LE_CORR_hybrid_Penman) -  rmse( LE_CORR,LE_CORR_RF),
							H2 =  rmse( LE_CORR,LE_CORR_hybrid_PT)-  rmse( LE_CORR,LE_CORR_RF),
							H3 =  rmse( LE_CORR,LE_CORR_hybrid_PMgs)-  rmse( LE_CORR,LE_CORR_RF),
							H4 =  rmse( LE_CORR,LE_CORR_hybrid_PMkc)-  rmse( LE_CORR,LE_CORR_RF),
							H5 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_drh)-  rmse( LE_CORR,LE_CORR_RF),
							H6 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_Frh)-  rmse( LE_CORR,LE_CORR_RF)
		)%>%
		gather(key = "model", value = "T_97")
)

#Fpar
aa <- aa %>% left_join(
	test_set2 %>% filter(Fpar < FPAR_03) %>%
		filter(site != "US-CS2") %>%
		summarise(H1 =  rmse( LE_CORR,LE_CORR_hybrid_Penman) -  rmse( LE_CORR,LE_CORR_RF),
							H2 =  rmse( LE_CORR,LE_CORR_hybrid_PT)-  rmse( LE_CORR,LE_CORR_RF),
							H3 =  rmse( LE_CORR,LE_CORR_hybrid_PMgs)-  rmse( LE_CORR,LE_CORR_RF),
							H4 =  rmse( LE_CORR,LE_CORR_hybrid_PMkc)-  rmse( LE_CORR,LE_CORR_RF),
							H5 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_drh)-  rmse( LE_CORR,LE_CORR_RF),
							H6 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_Frh)-  rmse( LE_CORR,LE_CORR_RF)
		)%>%
		gather(key = "model", value = "FPAR_03")
)

aa <- aa %>% left_join(
	test_set2 %>% filter(Fpar > FPAR_97) %>%
		filter(site != "US-CS2") %>%
		summarise(H1 =  rmse( LE_CORR,LE_CORR_hybrid_Penman) -  rmse( LE_CORR,LE_CORR_RF),
							H2 =  rmse( LE_CORR,LE_CORR_hybrid_PT)-  rmse( LE_CORR,LE_CORR_RF),
							H3 =  rmse( LE_CORR,LE_CORR_hybrid_PMgs)-  rmse( LE_CORR,LE_CORR_RF),
							H4 =  rmse( LE_CORR,LE_CORR_hybrid_PMkc)-  rmse( LE_CORR,LE_CORR_RF),
							H5 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_drh)-  rmse( LE_CORR,LE_CORR_RF),
							H6 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_Frh)-  rmse( LE_CORR,LE_CORR_RF)
		)%>%
		gather(key = "model", value = "FPAR_97")
)

#SW
aa <- aa %>% left_join(
	test_set2 %>% filter(SW_IN_F_MDS < SW_IN_03) %>%
		filter(site != "US-CS2") %>%
		summarise(H1 =  rmse( LE_CORR,LE_CORR_hybrid_Penman) -  rmse( LE_CORR,LE_CORR_RF),
							H2 =  rmse( LE_CORR,LE_CORR_hybrid_PT)-  rmse( LE_CORR,LE_CORR_RF),
							H3 =  rmse( LE_CORR,LE_CORR_hybrid_PMgs)-  rmse( LE_CORR,LE_CORR_RF),
							H4 =  rmse( LE_CORR,LE_CORR_hybrid_PMkc)-  rmse( LE_CORR,LE_CORR_RF),
							H5 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_drh)-  rmse( LE_CORR,LE_CORR_RF),
							H6 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_Frh)-  rmse( LE_CORR,LE_CORR_RF)
		)%>%
		gather(key = "model", value = "SW_IN_03")
)

aa <- aa %>% left_join(
	test_set2 %>% filter(SW_IN_F_MDS > SW_IN_97) %>%
		filter(site != "US-CS2") %>%
		summarise(H1 =  rmse( LE_CORR,LE_CORR_hybrid_Penman) -  rmse( LE_CORR,LE_CORR_RF),
							H2 =  rmse( LE_CORR,LE_CORR_hybrid_PT)-  rmse( LE_CORR,LE_CORR_RF),
							H3 =  rmse( LE_CORR,LE_CORR_hybrid_PMgs)-  rmse( LE_CORR,LE_CORR_RF),
							H4 =  rmse( LE_CORR,LE_CORR_hybrid_PMkc)-  rmse( LE_CORR,LE_CORR_RF),
							H5 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_drh)-  rmse( LE_CORR,LE_CORR_RF),
							H6 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_Frh)-  rmse( LE_CORR,LE_CORR_RF)
		)%>%
		gather(key = "model", value = "SW_IN_97")
)


#ga_H
aa <- aa %>% left_join(
	test_set2 %>% filter(ga_H < ga_H_03) %>%
		filter(site != "US-CS2") %>%
		summarise(H1 =  rmse( LE_CORR,LE_CORR_hybrid_Penman) -  rmse( LE_CORR,LE_CORR_RF),
							H2 =  rmse( LE_CORR,LE_CORR_hybrid_PT)-  rmse( LE_CORR,LE_CORR_RF),
							H3 =  rmse( LE_CORR,LE_CORR_hybrid_PMgs)-  rmse( LE_CORR,LE_CORR_RF),
							H4 =  rmse( LE_CORR,LE_CORR_hybrid_PMkc)-  rmse( LE_CORR,LE_CORR_RF),
							H5 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_drh)-  rmse( LE_CORR,LE_CORR_RF),
							H6 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_Frh)-  rmse( LE_CORR,LE_CORR_RF)
		)%>%
		gather(key = "model", value = "ga_H_03")
)

aa <- aa %>% left_join(
	test_set2 %>% filter(ga_H > ga_H_97) %>%
		filter(site != "US-CS2") %>%
		summarise(H1 =  rmse( LE_CORR,LE_CORR_hybrid_Penman) -  rmse( LE_CORR,LE_CORR_RF),
							H2 =  rmse( LE_CORR,LE_CORR_hybrid_PT)-  rmse( LE_CORR,LE_CORR_RF),
							H3 =  rmse( LE_CORR,LE_CORR_hybrid_PMgs)-  rmse( LE_CORR,LE_CORR_RF),
							H4 =  rmse( LE_CORR,LE_CORR_hybrid_PMkc)-  rmse( LE_CORR,LE_CORR_RF),
							H5 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_drh)-  rmse( LE_CORR,LE_CORR_RF),
							H6 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_Frh)-  rmse( LE_CORR,LE_CORR_RF)
		)%>%
		gather(key = "model", value = "ga_H_97")
)

#CO2
aa <- aa %>% left_join(
	test_set2 %>% filter(CO2_F_MDS < CO2_03) %>%
		filter(site != "US-CS2") %>%
		summarise(H1 =  rmse( LE_CORR,LE_CORR_hybrid_Penman) -  rmse( LE_CORR,LE_CORR_RF),
							H2 =  rmse( LE_CORR,LE_CORR_hybrid_PT)-  rmse( LE_CORR,LE_CORR_RF),
							H3 =  rmse( LE_CORR,LE_CORR_hybrid_PMgs)-  rmse( LE_CORR,LE_CORR_RF),
							H4 =  rmse( LE_CORR,LE_CORR_hybrid_PMkc)-  rmse( LE_CORR,LE_CORR_RF),
							H5 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_drh)-  rmse( LE_CORR,LE_CORR_RF),
							H6 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_Frh)-  rmse( LE_CORR,LE_CORR_RF)
		)%>%
		gather(key = "model", value = "CO2_03")
)

aa <- aa %>% left_join(
	test_set2 %>% filter(CO2_F_MDS > CO2_97) %>%
		filter(site != "US-CS2") %>%
		summarise(H1 =  rmse( LE_CORR,LE_CORR_hybrid_Penman) -  rmse( LE_CORR,LE_CORR_RF),
							H2 =  rmse( LE_CORR,LE_CORR_hybrid_PT)-  rmse( LE_CORR,LE_CORR_RF),
							H3 =  rmse( LE_CORR,LE_CORR_hybrid_PMgs)-  rmse( LE_CORR,LE_CORR_RF),
							H4 =  rmse( LE_CORR,LE_CORR_hybrid_PMkc)-  rmse( LE_CORR,LE_CORR_RF),
							H5 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_drh)-  rmse( LE_CORR,LE_CORR_RF),
							H6 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_Frh)-  rmse( LE_CORR,LE_CORR_RF)
		)%>%
		gather(key = "model", value = "CO2_97")
)


#VegWC
aa <- aa %>% left_join(
	test_set2 %>% filter(VegWC < VegWC_03) %>%
		filter(site != "US-CS2") %>%
		summarise(H1 =  rmse( LE_CORR,LE_CORR_hybrid_Penman) -  rmse( LE_CORR,LE_CORR_RF),
							H2 =  rmse( LE_CORR,LE_CORR_hybrid_PT)-  rmse( LE_CORR,LE_CORR_RF),
							H3 =  rmse( LE_CORR,LE_CORR_hybrid_PMgs)-  rmse( LE_CORR,LE_CORR_RF),
							H4 =  rmse( LE_CORR,LE_CORR_hybrid_PMkc)-  rmse( LE_CORR,LE_CORR_RF),
							H5 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_drh)-  rmse( LE_CORR,LE_CORR_RF),
							H6 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_Frh)-  rmse( LE_CORR,LE_CORR_RF)
		)%>%
		gather(key = "model", value = "VegWC_03")
)

aa <- aa %>% left_join(
	test_set2 %>% filter(VegWC > VegWC_97) %>%
		filter(site != "US-CS2") %>%
		summarise(H1 =  rmse( LE_CORR,LE_CORR_hybrid_Penman) -  rmse( LE_CORR,LE_CORR_RF),
							H2 =  rmse( LE_CORR,LE_CORR_hybrid_PT)-  rmse( LE_CORR,LE_CORR_RF),
							H3 =  rmse( LE_CORR,LE_CORR_hybrid_PMgs)-  rmse( LE_CORR,LE_CORR_RF),
							H4 =  rmse( LE_CORR,LE_CORR_hybrid_PMkc)-  rmse( LE_CORR,LE_CORR_RF),
							H5 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_drh)-  rmse( LE_CORR,LE_CORR_RF),
							H6 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_Frh)-  rmse( LE_CORR,LE_CORR_RF)
		)%>%
		gather(key = "model", value = "VegWC_97")
)



p1 <- aa %>%
	gather(key = "var",value = "rmse", -c("model")) %>%
	ggplot(aes(model,rmse)) +
	theme_bw()+
	geom_boxplot(outlier.alpha = 0) +
	geom_hline(yintercept = 0) +
	geom_jitter(aes(fill = var, shape = var), width=0.15, size = 2) +
	xlab("") + 
	ylab(expression("RMSE difference (Hybrid - pure ML) "*" (W "*m^-2*")")) +
	scale_shape_manual(name = "Extremes",
										 values = rep(c(25,24),8), 
										 labels = c("CO2 < 1%", "CO2 > 99%","FPAR < 1%", "FPAR > 99%", 
										 					 "gaH < 1%", "gaH > 99%","RH < 1%", "RH > 99%", 
										 					 "Rg < 1%", "Rg > 99%", "T < 1%", "T > 99%", 
										 					 "VegWC < 1%", "VegWC > 99%", "SM < 1%", "SM > 99%"),
										 guide = guide_legend(ncol = 2)
	) +
	scale_fill_manual(name = "Extremes",
										values = c("#999999","#999999","#E69F00","#E69F00", 
															 "#56B4E9","#56B4E9", "#009E73","#009E73", 
															 "#F0E442","#F0E442", "#D55E00","#D55E00", 
															 "#CC79A7","#CC79A7","#0072B2", "#0072B2"), 
										labels = c("CO2 < 1%", "CO2 > 99%","FPAR < 1%", "FPAR > 99%", 
															 "gaH < 1%", "gaH > 99%","RH < 1%", "RH > 99%", 
															 "Rg < 1%", "Rg > 99%", "T < 1%", "T > 99%", 
															 "VegWC < 1%", "VegWC > 99%", "SM < 1%", "SM > 99%"),
										guide = guide_legend(ncol = 2)
	)+
	scale_x_discrete(
		labels = c("H1" = "H-Penman", "H2" = "H-PT", "H3" = expression("H-PM-"*g[s]),
							 "H4" = expression("H-PM-"*K[c]), "H5" = expression("H-"*PM[RH]*"-"*d[RH]), "H6" = expression("H-"*PM[RH]*"-"*F[RH])
		)
	) +
	theme(axis.text.x = element_text(size =9))



#### 3%
extreme <- df %>% 
	group_by(site)%>%
	summarise(VWC_03 = quantile(soil_moisture,0.03),
						RH_03 = quantile(rha,0.03),
						T_03 = quantile(TA_F,0.03),
						FPAR_03 = quantile(Fpar,0.03),
						SW_IN_03 = quantile(SW_IN_F_MDS,0.03),
						ga_H_03 = quantile(ga_H,0.03),
						VegWC_03 = quantile(VegWC, 0.03),
						CO2_03 = quantile(CO2_F_MDS, 0.03),
						
						VWC_97 = quantile(soil_moisture,0.97),
						RH_97 = quantile(rha,0.97),
						T_97 = quantile(TA_F,0.97),
						FPAR_97 = quantile(Fpar,0.97),
						SW_IN_97 = quantile(SW_IN_F_MDS,0.97),
						ga_H_97 = quantile(ga_H,0.97),
						VegWC_97 = quantile(VegWC, 0.97),
						CO2_97 = quantile(CO2_F_MDS, 0.97)
	)	

df2 <- df %>% 
	left_join(extreme)

test_set2 <- df2[-index,]

# soil moisture
aa <- test_set2 %>% filter(soil_moisture < VWC_03) %>%
	filter(site != "US-CS2") %>%
	summarise(H1 =  rmse( LE_CORR,LE_CORR_hybrid_Penman) -  rmse( LE_CORR,LE_CORR_RF),
						H2 =  rmse( LE_CORR,LE_CORR_hybrid_PT)-  rmse( LE_CORR,LE_CORR_RF),
						H3 =  rmse( LE_CORR,LE_CORR_hybrid_PMgs)-  rmse( LE_CORR,LE_CORR_RF),
						H4 =  rmse( LE_CORR,LE_CORR_hybrid_PMkc)-  rmse( LE_CORR,LE_CORR_RF),
						H5 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_drh)-  rmse( LE_CORR,LE_CORR_RF),
						H6 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_Frh)-  rmse( LE_CORR,LE_CORR_RF)
	)%>%
	gather(key = "model", value = "VWC_03")

aa <- aa %>% left_join(
	test_set2 %>% filter(soil_moisture > VWC_97) %>%
		filter(site != "US-CS2") %>%
		summarise(H1 =  rmse( LE_CORR,LE_CORR_hybrid_Penman) -  rmse( LE_CORR,LE_CORR_RF),
							H2 =  rmse( LE_CORR,LE_CORR_hybrid_PT)-  rmse( LE_CORR,LE_CORR_RF),
							H3 =  rmse( LE_CORR,LE_CORR_hybrid_PMgs)-  rmse( LE_CORR,LE_CORR_RF),
							H4 =  rmse( LE_CORR,LE_CORR_hybrid_PMkc)-  rmse( LE_CORR,LE_CORR_RF),
							H5 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_drh)-  rmse( LE_CORR,LE_CORR_RF),
							H6 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_Frh)-  rmse( LE_CORR,LE_CORR_RF)
		)%>%
		gather(key = "model", value = "VWC_97")
)

#relative humidity
aa <- aa %>% left_join(
	test_set2 %>% filter(rha < RH_03) %>%
		filter(site != "US-CS2") %>%
		summarise(H1 =  rmse( LE_CORR,LE_CORR_hybrid_Penman) -  rmse( LE_CORR,LE_CORR_RF),
							H2 =  rmse( LE_CORR,LE_CORR_hybrid_PT)-  rmse( LE_CORR,LE_CORR_RF),
							H3 =  rmse( LE_CORR,LE_CORR_hybrid_PMgs)-  rmse( LE_CORR,LE_CORR_RF),
							H4 =  rmse( LE_CORR,LE_CORR_hybrid_PMkc)-  rmse( LE_CORR,LE_CORR_RF),
							H5 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_drh)-  rmse( LE_CORR,LE_CORR_RF),
							H6 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_Frh)-  rmse( LE_CORR,LE_CORR_RF)
		)%>%
		gather(key = "model", value = "RH_03")
)

aa <- aa %>% left_join(
	test_set2 %>% filter(rha > RH_97) %>%
		filter(site != "US-CS2") %>%
		summarise(H1 =  rmse( LE_CORR,LE_CORR_hybrid_Penman) -  rmse( LE_CORR,LE_CORR_RF),
							H2 =  rmse( LE_CORR,LE_CORR_hybrid_PT)-  rmse( LE_CORR,LE_CORR_RF),
							H3 =  rmse( LE_CORR,LE_CORR_hybrid_PMgs)-  rmse( LE_CORR,LE_CORR_RF),
							H4 =  rmse( LE_CORR,LE_CORR_hybrid_PMkc)-  rmse( LE_CORR,LE_CORR_RF),
							H5 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_drh)-  rmse( LE_CORR,LE_CORR_RF),
							H6 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_Frh)-  rmse( LE_CORR,LE_CORR_RF)
		)%>%
		gather(key = "model", value = "RH_97")
)

#temperature
aa <- aa %>% left_join(
	test_set2 %>% filter(TA_F < T_03) %>%
		filter(site != "US-CS2") %>%
		summarise(H1 =  rmse( LE_CORR,LE_CORR_hybrid_Penman) -  rmse( LE_CORR,LE_CORR_RF),
							H2 =  rmse( LE_CORR,LE_CORR_hybrid_PT)-  rmse( LE_CORR,LE_CORR_RF),
							H3 =  rmse( LE_CORR,LE_CORR_hybrid_PMgs)-  rmse( LE_CORR,LE_CORR_RF),
							H4 =  rmse( LE_CORR,LE_CORR_hybrid_PMkc)-  rmse( LE_CORR,LE_CORR_RF),
							H5 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_drh)-  rmse( LE_CORR,LE_CORR_RF),
							H6 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_Frh)-  rmse( LE_CORR,LE_CORR_RF)
		)%>%
		gather(key = "model", value = "T_03")
)

aa <- aa %>% left_join(
	test_set2 %>% filter(TA_F > T_97) %>%
		filter(site != "US-CS2") %>%
		summarise(H1 =  rmse( LE_CORR,LE_CORR_hybrid_Penman) -  rmse( LE_CORR,LE_CORR_RF),
							H2 =  rmse( LE_CORR,LE_CORR_hybrid_PT)-  rmse( LE_CORR,LE_CORR_RF),
							H3 =  rmse( LE_CORR,LE_CORR_hybrid_PMgs)-  rmse( LE_CORR,LE_CORR_RF),
							H4 =  rmse( LE_CORR,LE_CORR_hybrid_PMkc)-  rmse( LE_CORR,LE_CORR_RF),
							H5 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_drh)-  rmse( LE_CORR,LE_CORR_RF),
							H6 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_Frh)-  rmse( LE_CORR,LE_CORR_RF)
		)%>%
		gather(key = "model", value = "T_97")
)

#Fpar
aa <- aa %>% left_join(
	test_set2 %>% filter(Fpar < FPAR_03) %>%
		filter(site != "US-CS2") %>%
		summarise(H1 =  rmse( LE_CORR,LE_CORR_hybrid_Penman) -  rmse( LE_CORR,LE_CORR_RF),
							H2 =  rmse( LE_CORR,LE_CORR_hybrid_PT)-  rmse( LE_CORR,LE_CORR_RF),
							H3 =  rmse( LE_CORR,LE_CORR_hybrid_PMgs)-  rmse( LE_CORR,LE_CORR_RF),
							H4 =  rmse( LE_CORR,LE_CORR_hybrid_PMkc)-  rmse( LE_CORR,LE_CORR_RF),
							H5 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_drh)-  rmse( LE_CORR,LE_CORR_RF),
							H6 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_Frh)-  rmse( LE_CORR,LE_CORR_RF)
		)%>%
		gather(key = "model", value = "FPAR_03")
)

aa <- aa %>% left_join(
	test_set2 %>% filter(Fpar > FPAR_97) %>%
		filter(site != "US-CS2") %>%
		summarise(H1 =  rmse( LE_CORR,LE_CORR_hybrid_Penman) -  rmse( LE_CORR,LE_CORR_RF),
							H2 =  rmse( LE_CORR,LE_CORR_hybrid_PT)-  rmse( LE_CORR,LE_CORR_RF),
							H3 =  rmse( LE_CORR,LE_CORR_hybrid_PMgs)-  rmse( LE_CORR,LE_CORR_RF),
							H4 =  rmse( LE_CORR,LE_CORR_hybrid_PMkc)-  rmse( LE_CORR,LE_CORR_RF),
							H5 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_drh)-  rmse( LE_CORR,LE_CORR_RF),
							H6 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_Frh)-  rmse( LE_CORR,LE_CORR_RF)
		)%>%
		gather(key = "model", value = "FPAR_97")
)

#SW
aa <- aa %>% left_join(
	test_set2 %>% filter(SW_IN_F_MDS < SW_IN_03) %>%
		filter(site != "US-CS2") %>%
		summarise(H1 =  rmse( LE_CORR,LE_CORR_hybrid_Penman) -  rmse( LE_CORR,LE_CORR_RF),
							H2 =  rmse( LE_CORR,LE_CORR_hybrid_PT)-  rmse( LE_CORR,LE_CORR_RF),
							H3 =  rmse( LE_CORR,LE_CORR_hybrid_PMgs)-  rmse( LE_CORR,LE_CORR_RF),
							H4 =  rmse( LE_CORR,LE_CORR_hybrid_PMkc)-  rmse( LE_CORR,LE_CORR_RF),
							H5 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_drh)-  rmse( LE_CORR,LE_CORR_RF),
							H6 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_Frh)-  rmse( LE_CORR,LE_CORR_RF)
		)%>%
		gather(key = "model", value = "SW_IN_03")
)

aa <- aa %>% left_join(
	test_set2 %>% filter(SW_IN_F_MDS > SW_IN_97) %>%
		filter(site != "US-CS2") %>%
		summarise(H1 =  rmse( LE_CORR,LE_CORR_hybrid_Penman) -  rmse( LE_CORR,LE_CORR_RF),
							H2 =  rmse( LE_CORR,LE_CORR_hybrid_PT)-  rmse( LE_CORR,LE_CORR_RF),
							H3 =  rmse( LE_CORR,LE_CORR_hybrid_PMgs)-  rmse( LE_CORR,LE_CORR_RF),
							H4 =  rmse( LE_CORR,LE_CORR_hybrid_PMkc)-  rmse( LE_CORR,LE_CORR_RF),
							H5 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_drh)-  rmse( LE_CORR,LE_CORR_RF),
							H6 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_Frh)-  rmse( LE_CORR,LE_CORR_RF)
		)%>%
		gather(key = "model", value = "SW_IN_97")
)


#ga_H
aa <- aa %>% left_join(
	test_set2 %>% filter(ga_H < ga_H_03) %>%
		filter(site != "US-CS2") %>%
		summarise(H1 =  rmse( LE_CORR,LE_CORR_hybrid_Penman) -  rmse( LE_CORR,LE_CORR_RF),
							H2 =  rmse( LE_CORR,LE_CORR_hybrid_PT)-  rmse( LE_CORR,LE_CORR_RF),
							H3 =  rmse( LE_CORR,LE_CORR_hybrid_PMgs)-  rmse( LE_CORR,LE_CORR_RF),
							H4 =  rmse( LE_CORR,LE_CORR_hybrid_PMkc)-  rmse( LE_CORR,LE_CORR_RF),
							H5 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_drh)-  rmse( LE_CORR,LE_CORR_RF),
							H6 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_Frh)-  rmse( LE_CORR,LE_CORR_RF)
		)%>%
		gather(key = "model", value = "ga_H_03")
)

aa <- aa %>% left_join(
	test_set2 %>% filter(ga_H > ga_H_97) %>%
		filter(site != "US-CS2") %>%
		summarise(H1 =  rmse( LE_CORR,LE_CORR_hybrid_Penman) -  rmse( LE_CORR,LE_CORR_RF),
							H2 =  rmse( LE_CORR,LE_CORR_hybrid_PT)-  rmse( LE_CORR,LE_CORR_RF),
							H3 =  rmse( LE_CORR,LE_CORR_hybrid_PMgs)-  rmse( LE_CORR,LE_CORR_RF),
							H4 =  rmse( LE_CORR,LE_CORR_hybrid_PMkc)-  rmse( LE_CORR,LE_CORR_RF),
							H5 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_drh)-  rmse( LE_CORR,LE_CORR_RF),
							H6 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_Frh)-  rmse( LE_CORR,LE_CORR_RF)
		)%>%
		gather(key = "model", value = "ga_H_97")
)

#CO2
aa <- aa %>% left_join(
	test_set2 %>% filter(CO2_F_MDS < CO2_03) %>%
		filter(site != "US-CS2") %>%
		summarise(H1 =  rmse( LE_CORR,LE_CORR_hybrid_Penman) -  rmse( LE_CORR,LE_CORR_RF),
							H2 =  rmse( LE_CORR,LE_CORR_hybrid_PT)-  rmse( LE_CORR,LE_CORR_RF),
							H3 =  rmse( LE_CORR,LE_CORR_hybrid_PMgs)-  rmse( LE_CORR,LE_CORR_RF),
							H4 =  rmse( LE_CORR,LE_CORR_hybrid_PMkc)-  rmse( LE_CORR,LE_CORR_RF),
							H5 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_drh)-  rmse( LE_CORR,LE_CORR_RF),
							H6 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_Frh)-  rmse( LE_CORR,LE_CORR_RF)
		)%>%
		gather(key = "model", value = "CO2_03")
)

aa <- aa %>% left_join(
	test_set2 %>% filter(CO2_F_MDS > CO2_97) %>%
		filter(site != "US-CS2") %>%
		summarise(H1 =  rmse( LE_CORR,LE_CORR_hybrid_Penman) -  rmse( LE_CORR,LE_CORR_RF),
							H2 =  rmse( LE_CORR,LE_CORR_hybrid_PT)-  rmse( LE_CORR,LE_CORR_RF),
							H3 =  rmse( LE_CORR,LE_CORR_hybrid_PMgs)-  rmse( LE_CORR,LE_CORR_RF),
							H4 =  rmse( LE_CORR,LE_CORR_hybrid_PMkc)-  rmse( LE_CORR,LE_CORR_RF),
							H5 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_drh)-  rmse( LE_CORR,LE_CORR_RF),
							H6 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_Frh)-  rmse( LE_CORR,LE_CORR_RF)
		)%>%
		gather(key = "model", value = "CO2_97")
)


#VegWC
aa <- aa %>% left_join(
	test_set2 %>% filter(VegWC < VegWC_03) %>%
		filter(site != "US-CS2") %>%
		summarise(H1 =  rmse( LE_CORR,LE_CORR_hybrid_Penman) -  rmse( LE_CORR,LE_CORR_RF),
							H2 =  rmse( LE_CORR,LE_CORR_hybrid_PT)-  rmse( LE_CORR,LE_CORR_RF),
							H3 =  rmse( LE_CORR,LE_CORR_hybrid_PMgs)-  rmse( LE_CORR,LE_CORR_RF),
							H4 =  rmse( LE_CORR,LE_CORR_hybrid_PMkc)-  rmse( LE_CORR,LE_CORR_RF),
							H5 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_drh)-  rmse( LE_CORR,LE_CORR_RF),
							H6 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_Frh)-  rmse( LE_CORR,LE_CORR_RF)
		)%>%
		gather(key = "model", value = "VegWC_03")
)

aa <- aa %>% left_join(
	test_set2 %>% filter(VegWC > VegWC_97) %>%
		filter(site != "US-CS2") %>%
		summarise(H1 =  rmse( LE_CORR,LE_CORR_hybrid_Penman) -  rmse( LE_CORR,LE_CORR_RF),
							H2 =  rmse( LE_CORR,LE_CORR_hybrid_PT)-  rmse( LE_CORR,LE_CORR_RF),
							H3 =  rmse( LE_CORR,LE_CORR_hybrid_PMgs)-  rmse( LE_CORR,LE_CORR_RF),
							H4 =  rmse( LE_CORR,LE_CORR_hybrid_PMkc)-  rmse( LE_CORR,LE_CORR_RF),
							H5 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_drh)-  rmse( LE_CORR,LE_CORR_RF),
							H6 =  rmse( LE_CORR,LE_CORR_hybrid_PMRH_Frh)-  rmse( LE_CORR,LE_CORR_RF)
		)%>%
		gather(key = "model", value = "VegWC_97")
)



p2 <- aa %>%
	gather(key = "var",value = "rmse", -c("model")) %>%
	ggplot(aes(model,rmse)) +
	theme_bw()+
	geom_boxplot(outlier.alpha = 0) +
	geom_hline(yintercept = 0) +
	geom_jitter(aes(fill = var, shape = var), width=0.15, size = 2) +
	xlab("") + 
	ylab(expression("RMSE difference (Hybrid - pure ML) "*" (W "*m^-2*")")) +
	scale_shape_manual(name = "Extremes",
										 values = rep(c(25,24),8), 
										 labels = c("CO2 < 3%", "CO2 > 97%","FPAR < 3%", "FPAR > 97%", 
										 					 "gaH < 3%", "gaH > 97%","RH < 3%", "RH > 97%", 
										 					 "Rg < 3%", "Rg > 97%", "T < 3%", "T > 97%", 
										 					 "VegWC < 3%", "VegWC > 97%", "SM < 3%", "SM > 97%"),
										 guide = guide_legend(ncol = 2)
	) +
	scale_fill_manual(name = "Extremes",
										values = c("#999999","#999999","#E69F00","#E69F00", 
															 "#56B4E9","#56B4E9", "#009E73","#009E73", 
															 "#F0E442","#F0E442", "#D55E00","#D55E00", 
															 "#CC79A7","#CC79A7","#0072B2", "#0072B2"), 
										labels = c("CO2 < 3%", "CO2 > 97%","FPAR < 3%", "FPAR > 97%", 
															 "gaH < 3%", "gaH > 97%","RH < 3%", "RH > 97%", 
															 "Rg < 3%", "Rg > 97%", "T < 3%", "T > 97%", 
															 "VegWC < 3%", "VegWC > 97%", "SM < 3%", "SM > 97%"),
										guide = guide_legend(ncol = 2)
	)+
	scale_x_discrete(
		labels = c("H1" = "H-Penman", "H2" = "H-PT", "H3" = expression("H-PM-"*g[s]),
							 "H4" = expression("H-PM-"*K[c]), "H5" = expression("H-"*PM[RH]*"-"*d[RH]), "H6" = expression("H-"*PM[RH]*"-"*F[RH])
		)
	) +
	theme(axis.text.x = element_text(size =9))


tiff(filename = "./Figures/tif_Figure7.tif",width = 7.5, height = 4,units = "in", res = 300)
p2
dev.off()



#####################################
######### variable importance
#####################################
# direct LE
imp  <- varImp(RF_LE_CORR, scale = F)$importance
imp$varnames <- c("FPAR","SM","VegWC","T","CO2","Rg","gaH","RH","AE")
imp$category <- c("Satellite","Satellite","Satellite","Meteo","Meteo","Meteo","Meteo","Meteo","Meteo")
a <- imp %>% ggplot(aes(x=reorder(varnames, Overall),fill=category, weight=Overall )) +
	geom_bar(width = 0.6) +
	theme_bw() +
	ylab("Importance (%IncMSE)") +
	theme(axis.title.y = element_blank(),
				axis.title.x = element_text(size = 7),
				title = element_text(size = 8),
				axis.text = element_text(size = 7)
				) +
	coord_flip() + scale_fill_manual(values = c("#0072B2","#D55E00")) +
	ggtitle(expression(italic("Pure-ML")))

# Penman
imp  <- varImp(RF_Penman_stress_cor, scale = F)$importance
imp$varnames <- c("FPAR","SM","VegWC","T","CO2","Rg","gaH","RH")
imp$category <- c("Satellite","Satellite","Satellite","Meteo","Meteo","Meteo","Meteo","Meteo")
b <- imp %>% ggplot(aes(x=reorder(varnames, Overall),fill=category, weight=Overall )) +
	geom_bar(width = 0.6) +
	theme_bw() +
	ylab("Importance (%IncMSE)") +
	theme(axis.title.y = element_blank(),
				axis.title.x = element_text(size = 7),
				title = element_text(size = 8),
				axis.text = element_text(size = 7)
	) +
	coord_flip() + scale_fill_manual(values = c("#0072B2","#D55E00")) +
	ggtitle(expression(italic("H-Penman")))

# PT
imp  <- varImp(RF_PT_stress_cor, scale = F)$importance
imp$varnames <- c("FPAR","SM","VegWC","T","CO2","Rg","gaH","RH")
imp$category <- c("Satellite","Satellite","Satellite","Meteo","Meteo","Meteo","Meteo","Meteo")
c <- imp %>% ggplot(aes(x=reorder(varnames, Overall),fill=category, weight=Overall )) +
	geom_bar(width = 0.6) +
	theme_bw() +
	ylab("Importance (%IncMSE)") +
	theme(axis.title.y = element_blank(),
				axis.title.x = element_text(size = 7),
				title = element_text(size = 8),
				axis.text = element_text(size = 7)
	) +
	coord_flip() + scale_fill_manual(values = c("#0072B2","#D55E00")) +
	ggtitle(expression(italic("H-PT")))

# PM-gs
imp  <- varImp(RF_Gs_log_cor, scale = F)$importance
imp$varnames <- c("FPAR","SM","VegWC","T","CO2","Rg","gaH","RH")
imp$category <- c("Satellite","Satellite","Satellite","Meteo","Meteo","Meteo","Meteo","Meteo")
d <- imp %>% ggplot(aes(x=reorder(varnames, Overall),fill=category, weight=Overall )) +
	geom_bar(width = 0.6) +
	theme_bw() +
	ylab("Importance (%IncMSE)") +
	theme(axis.title.y = element_blank(),
				axis.title.x = element_text(size = 7),
				title = element_text(size = 8),
				axis.text = element_text(size = 7)
	) +
	coord_flip() + scale_fill_manual(values = c("#0072B2","#D55E00")) +
	ggtitle(expression(italic("H-PM-"*g[s])))

# PM-Kc
imp  <- varImp(RF_PM_stress_cor, scale = F)$importance
imp$varnames <- c("FPAR","SM","VegWC","T","CO2","Rg","gaH","RH")
imp$category <- c("Satellite","Satellite","Satellite","Meteo","Meteo","Meteo","Meteo","Meteo")
e <- imp %>% ggplot(aes(x=reorder(varnames, Overall),fill=category, weight=Overall )) +
	geom_bar(width = 0.6) +
	theme_bw() +
	ylab("Importance (%IncMSE)") +
	theme(axis.title.y = element_blank(),
				axis.title.x = element_text(size = 7),
				title = element_text(size = 8),
				axis.text = element_text(size = 7)
	) +
	coord_flip() + scale_fill_manual(values = c("#0072B2","#D55E00")) +
	ggtitle(expression(italic("H-PM-"*K[c])))

# PMRH-dRH
imp  <- varImp(RF_drh_cor, scale = F)$importance
imp$varnames <- c("FPAR","SM","VegWC","T","CO2","Rg","gaH","RH")
imp$category <- c("Satellite","Satellite","Satellite","Meteo","Meteo","Meteo","Meteo","Meteo")
f <- imp %>% ggplot(aes(x=reorder(varnames, Overall),fill=category, weight=Overall )) +
	geom_bar(width = 0.6) +
	theme_bw() +
	ylab("Importance (%IncMSE)") +
	theme(axis.title.y = element_blank(),
				axis.title.x = element_text(size = 7),
				title = element_text(size = 8),
				axis.text = element_text(size = 7)
	) +
	coord_flip() + scale_fill_manual(values = c("#0072B2","#D55E00")) +
	ggtitle(expression(italic("H-"*PM[RH]*"-"*d[RH])))

# PMRH-FRH
imp  <- varImp(RF_Frh_cor, scale = F)$importance
imp$varnames <- c("FPAR","SM","VegWC","T","CO2","Rg","gaH","RH")
imp$category <- c("Satellite","Satellite","Satellite","Meteo","Meteo","Meteo","Meteo","Meteo")
g <- imp %>% ggplot(aes(x=reorder(varnames, Overall),fill=category, weight=Overall )) +
	geom_bar(width = 0.6) +
	theme_bw() +
	ylab("Importance (%IncMSE)") +
	theme(axis.title.y = element_blank(),
				axis.title.x = element_text(size = 7),
				title = element_text(size = 8),
				axis.text = element_text(size = 7)
	) +
	coord_flip() + scale_fill_manual(values = c("#0072B2","#D55E00")) +
	ggtitle(expression(italic("H-"*PM[RH]*"-"*F[RH])))

p <- plot_grid(
  plot_grid(a + theme(legend.position = NaN),
            b + theme(legend.position = NaN),
            c + theme(legend.position = NaN),
            d + theme(legend.position = NaN), ncol = 4, labels = c("a","b","c","d")),
  plot_grid(NULL,e + theme(legend.position = NaN),
            f + theme(legend.position = NaN),
            g + theme(legend.position = NaN),
            get_legend(a),ncol = 5, rel_widths = c(0.5,1,1,1,0.5), labels = c("","e","f","g","")),
  ncol = 1 
)

tiff(filename = "./Figures/tif_Figure8.tif",width = 7.5, height = 5,units = "in", res = 300)
p
dev.off()

