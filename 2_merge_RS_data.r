rm(list = ls())

library(tidyverse)
library(lubridate)

# load data
df <- read_csv("./FLX_AA_DD_merged.csv")
df <- df %>% mutate(DATE = date(PDATE),
										YEAR = year(DATE)
										)
head(df)
info <- df %>% group_by(site) %>% tally()

###############################################################################
### 1. merge satellite data into EC data
###############################################################################

# satellite data (retrieve from google earth engine)
VWC <- read_csv("./SMAP_soil_moisture_data.csv")
VegWC <- read_csv("./SMAP_vegetation_water_content_data.csv")
Fpar <- read_csv("./MODIS_Fpar_data.csv")
# Fpar scaling
Fpar$Fpar <- Fpar$Fpar*0.01
head(Fpar)

# wrangling
VWC <- VWC %>% 
	mutate(DATE = date(date),
				 site = name
	) %>%
	select(!c(date,name))

VegWC <- VegWC %>% 
	mutate(DATE = date(date),
				 site = name
	) %>%
	select(!c(date,name))

Fpar <- Fpar %>% 
	mutate(DATE = date(date),
				 site = name
	) %>%
	select(!c(date,name))

### merge
df2 <- left_join(df,VWC)
df2 <- left_join(df2,VegWC)
df2 <- left_join(df2,Fpar)

### RS gap-fill (upto 5 days)
df2 <- df2 %>%
	group_by(site) %>%
	mutate(soil_moisture_raw = soil_moisture,
				 soil_moisture = ifelse(is.na(soil_moisture), 
															lead(soil_moisture, default = last(soil_moisture)), 
															soil_moisture),
				 soil_moisture = ifelse(is.na(soil_moisture), 
				 											 lead(soil_moisture, default = last(soil_moisture)), 
				 											 soil_moisture),
				 soil_moisture = ifelse(is.na(soil_moisture), 
				 											 lead(soil_moisture, default = last(soil_moisture)), 
				 											 soil_moisture),
				 soil_moisture = ifelse(is.na(soil_moisture), 
				 											 lead(soil_moisture, default = last(soil_moisture)), 
				 											 soil_moisture),
				 VegWC_raw = VegWC,
				 VegWC = ifelse(is.na(VegWC),
				 											 lead(VegWC, default = last(VegWC)),
				 											 VegWC),
				 VegWC = ifelse(is.na(VegWC),
				 											 lead(VegWC, default = last(VegWC)),
				 											 VegWC),
				 VegWC = ifelse(is.na(VegWC),
				 											 lead(VegWC, default = last(VegWC)),
				 											 VegWC),
				 VegWC = ifelse(is.na(VegWC),
				 							 lead(VegWC, default = last(VegWC)),
				 							 VegWC),
				 Fpar_raw = Fpar,
				 Fpar = ifelse(is.na(Fpar), 
				 											 lead(Fpar, default = last(Fpar)), 
				 											 Fpar),
				 Fpar = ifelse(is.na(Fpar), 
				 											 lead(Fpar, default = last(Fpar)), 
				 											 Fpar),
				 Fpar = ifelse(is.na(Fpar), 
				 											 lead(Fpar, default = last(Fpar)), 
				 											 Fpar),
				 Fpar = ifelse(is.na(Fpar), 
				 							lead(Fpar, default = last(Fpar)), 
				 							Fpar)
				 )


###############################################################################
### 2. EC data qc
###############################################################################

# quality flag
df2 <- df2 %>% filter(LE_F_MDS_QC >= 0.8) 
## variables which should not be NA
vars2 <- c("site","PDATE","TA_F", "rha", "PA_F", "CO2_F_MDS", "NETRAD","G_F_MDS","SW_IN_F_MDS",
					 "H_F_MDS", "H_CORR", "LE_F_MDS",	"LE_CORR","WS_F", "USTAR", "soil_moisture", "Fpar", "LC_Type1")
df2 <- df2 %>% drop_na(vars2)
# relative humidity qc
df2 <- df2 %>% filter(rha > 0 )

# remove unrealistic surface energy balance
df2 <- df2 %>% filter(abs(NETRAD - G_F_MDS - LE_CORR - H_CORR) < 100 )
df2 <- df2 %>% filter(NETRAD - G_F_MDS < 280 )
df2 <- df2 %>% filter(NETRAD - G_F_MDS > -50 )
# plot(df2$NETRAD - df2$G_F_MDS, df2$LE_CORR + df2$H_CORR)
# abline(a=0,b=1)

###############################
## save data
###############################
getwd()
write.csv(df, "FLX_RS_merged_filtered2.csv",row.names = F)

