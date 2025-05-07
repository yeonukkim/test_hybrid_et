rm(list = ls())

library("bigleaf")
library("tidyverse")
library("lubridate")

###############################################################################
####### 1. Merge daily EC data 
###############################################################################

## path
setwd("path to ameriflux data")
ls <- list.files()
DD_ls <- grep('FULLSET_DD', ls, value=TRUE)

## selected variables
vars <- c("site","PDATE","TA_F", "TA_F_QC", "VPD_F", "VPD_F_QC", "PA_F",	"P_F", "CO2_F_MDS", "CO2_F_MDS_QC",
					"NETRAD",	"NETRAD_QC", "G_F_MDS","G_F_MDS_QC",	"SW_IN_F_MDS","SW_IN_F_MDS_QC",
					"H_F_MDS", "H_F_MDS_QC", "H_CORR", "LE_F_MDS",	"LE_F_MDS_QC", "LE_CORR",
					"WS_F", "WS_F_QC", "USTAR", "USTAR_QC")

## test i == 1
df <- read.csv(DD_ls[1])
df[df==-9999] <- NA
# get date
YEAR <- floor(df$TIMESTAMP/1e4)
MONTH <- floor(df$TIMESTAMP/1e2)-YEAR*100
DATE <- df$TIMESTAMP%%100
PDATE <- ISOdate(YEAR,MONTH,DATE)
# get site name
site <- substring(DD_ls[1],5,10)
# define a new data frame 
LEdata <- cbind(site,PDATE,df)
LEdata <- LEdata %>% 
  #required variables
	dplyr::select(vars) %>%
  # filter data range based on SMAP availability
	filter(PDATE >= ISOdate(2015,4,1))
# resulting data frame
mergeddata <- LEdata

# loop to get daily data from every sites
for (i in 2: length(DD_ls)){
	df <- read.csv(DD_ls[i])
	df[df==-9999] <- NA

	#get date
	YEAR <- floor(df$TIMESTAMP/1e4)
	MONTH <- floor(df$TIMESTAMP/1e2)-YEAR*100
	DATE <- df$TIMESTAMP%%100
	PDATE <- ISOdate(YEAR,MONTH,DATE)
	#get site name
	site <- substring(DD_ls[i],5,10)
	# define a new data frame 
	LEdata <- cbind(site,PDATE,df)
	# filter data range based on SMAP availability
	LEdata <- LEdata %>% 
		filter(PDATE >= ISOdate(2015,4,1))
	
	# get column names
	colnames_temp <- colnames(LEdata)
	
	# add site only if all required variables are available
	if(!is.na(sum(match(vars, colnames_temp)))){
		LEdata <- LEdata %>% 
			dplyr::select(vars)
		# merge data
		mergeddata <- full_join(mergeddata,LEdata)
		}
}

###############################################################################
####### 2. Calculate additional variables for hybrid models
###############################################################################

df <- mergeddata

# available energy
df$AE <- df$NETRAD - df$G_F_MDS
# relative humidity
df$rha <- VPD.to.rH(df$VPD_F/10,df$TA_F, Esat.formula = "Allen_1998")
# specific heat capacity
cp <- bigleaf.constants()$cp
# air density
density <- air.density(df$TA_F,df$PA_F)
df$density <- density
# saturation vapor pressure
esat_air <- Esat.slope(df$TA_F)$Esat
df$esat_air <- esat_air
# specific humidity
q <- VPD.to.q(df$VPD_F/10,df$TA_F, df$PA_F, Esat.formula = "Allen_1998")
df$q <- q
# saturation vapor pressure slope
Delta <- Esat.slope(df$TA_F)$Delta
df$Delta <- Delta
# psychrometric constant
gamma <- psychrometric.constant(df$TA_F,df$PA_F)
df$gamma <- gamma
# latent heat of vaporization
Lv <- latent.heat.vaporization(df$TA_F)
df$Lv <- Lv

# aerodynamic resistance (conductance)
df$WS_F <- ifelse(df$WS_F < 0.2, 0.2, df$WS_F)
ra <- df$WS_F/(df$USTAR^2)
rb <- 6.2*df$USTAR^(-0.667)
ra_H <- ra + rb
ra_H <- ifelse(ra_H > 1000, 1000,ifelse(ra_H < 10, 10, ra_H))
df$ra_H <- ra_H
df$ga_H <- 1/ra_H

###############################################################################
####### 3. Calculate intermediate parameters (target of ML models)
###############################################################################

# method 1: Penman PE * stress
df$PET_Penman <- (Delta*df$AE+density*cp*df$VPD_F/10/ra_H)/(Delta+gamma)
df$Penman_stress_cor_raw <- ifelse(df$LE_CORR/df$PET_Penman < 0, 0, df$LE_CORR/df$PET_Penman)

# method 2: PT PE * stress
df$PET_PT <- 1.26*Delta/(Delta+gamma)*df$AE
df$PT_stress_cor_raw <- ifelse(df$LE_CORR/df$PET_PT < 0, 0, df$LE_CORR/df$PET_PT)

# method 3: surface conductance
Gs_cor <- df$LE_CORR*Ga_H*gamma/(Delta*(df$NETRAD - df$G_F_MDS) + density*cp*Ga_H*df$VPD_F/10 - df$LE_CORR*(Delta+gamma))
df <- df %>% mutate(Gs_cor_raw = ifelse(Gs_cor > 0, Gs_cor,NA))

# method 4: PM PE * Kc
df$PET_PM <- (Delta*(df$NETRAD - df$G_F_MDS)+density*cp*df$VPD_F/10/ra_H)/(Delta+gamma+gamma*70/ra_H)
df$PM_stress_cor_raw <- ifelse(df$LE_CORR/df$PET_PM < 0, 0, df$LE_CORR/df$PET_PM)

# method 5: Frh
df$Frh_cor <- (df$LE_CORR*(1+df$rha*Delta/gamma) - df$rha*Delta/gamma*(df$NETRAD - df$G_F_MDS))/(density*cp*esat_air)*gamma

# method 6: drh
df$drh_cor <- df$Frh_cor*ra_H

head(df)

##################################
## save data
##################################
getwd()
setwd("output path")
write.csv(df, "FLX_AA_DD_merged.csv",row.names = F)
