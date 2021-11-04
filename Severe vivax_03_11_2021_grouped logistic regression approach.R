#=================================================================
# Title		: Literature review of severe vivax malaria 
# Data version	: 11-Oct-2021 (from APP)
# Script Date	: 03-Nov-2021
#=================================================================
rm(list=ls())
library(readxl)

# meta-analysis packages
library(meta)
library(metafor)
library(dplyr)
library(lme4)

setwd("D:/_Colleagues/Aung Pyae Phyo/Rebuttal 2")

#=================================
# Read data 
#=================================
dat0 <- read_excel("S4 Table. List of articles included in meta-analyses_11Oct2021.xlsx", 
		sheet = "List of articles included"
		)

dat <- dat0[which(dat0$meta_analysis %in% c("Yes")),]

# number of unique articles
nrow(unique<-dat [which(!duplicated(dat$unique_ID)),])

# Add the number of cases based on WHO and other definitions to obtain overall number of cases
dat$severe_pv_all 	<- dat$Severe_WHO+dat$Severe_other
dat$cerebral_all 	<- dat$cerebral_WHO+dat$cerebral_other
dat$renal_all 		<- dat$renal_WHO+dat$renal_other
dat$respiratory_all 	<- dat$respiratory_WHO+dat$respiratory_other

#===========================================================
# Exclude studies in pregnancies for the main meta-analysis
#===========================================================
# Pregnancy studies only
preg <- dat[which(dat$pregnancy=="Yes"),] 

# Exclude pregnancy studies for main meta-analysis
dat <- dat[which(dat$pregnancy!="Yes"),]
dat <- droplevels(dat)

# number of unique articles
nrow(unique <- dat[which(!duplicated(dat$unique_ID)),])

#=====================================================================================================
# Exclude studies with either missing numerator or denominator for estimating severe vivax proportion
#=====================================================================================================
dat1 <- dat[which(!is.na(dat$severe_pv_all) & !is.na(dat$sample_size)),]

dat1 %>% 
	dplyr::summarise(
		n_pub= length(unique(unique_ID)),
		n_severe = sum(severe_pv_all, na.rm=T),
		n_total= sum(sample_size, na.rm=T),
		paste = paste(sum(severe_pv_all, na.rm=T), sum(sample_size, na.rm=T), sep="/")
	)

dat1 %>% 
	dplyr::group_by(region) %>%
	dplyr::summarise(
		n_pub= length(unique(unique_ID)),
		n_severe = sum(severe_pv_all, na.rm=T),
		n_total= sum(sample_size, na.rm=T),
		paste = paste(sum(severe_pv_all, na.rm=T), sum(sample_size, na.rm=T), sep="/")
	)

dat1 %>% 
	dplyr::group_by(settings) %>%
	dplyr::summarise(
		n_pub= length(unique(unique_ID)),
		n_severe = sum(severe_pv_all, na.rm=T),
		n_total= sum(sample_size, na.rm=T),
		paste = paste(sum(severe_pv_all, na.rm=T), sum(sample_size, na.rm=T), sep="/")
	)

#################################################################################
# Meta-analysis I: Estimating proportion of patients with severe vivax using ALL
#################################################################################
# Meta-analysis of proportion
(meta.prop <- metaprop(
			data = dat1,
			severe_pv_all,
			sample_size, 
			prediction=TRUE,
			studlab = unique_ID			)
		)
#========================================
# Subgroup analysis by region & settings
#========================================
update.meta(meta.prop,byvar=region_1)
update.meta(meta.prop,byvar=settings_1)


#=========================================================================
# Mixed effects logistic regression appraoch: grouped logistic regression
#=========================================================================
table(dat1$region,dat1$settings)

# Combine the levels of region and hopsitalisation to avoid rank deficiency in fixed effects matrix when fitting
# interaction between region and settings

dat1$region_1 <- dat1$region
dat1$region_1 [dat1$region_1=="Oceania"] <- "Asia"

dat1$settings_1 <- dat1$settings
dat1$settings_1 [dat1$settings_1=="Other"] <- "outpatients"

# new matrix
table(dat1$region_1,dat1$settings_1)

#=========================================================================
# Mixed effects logistic regression appraoch: ALL DEFINITION
#=========================================================================

#------------------
# Baseline model
#------------------
model_0 <- glmer(cbind(severe_pv_all,sample_size - severe_pv_all) ~    (1 | unique_ID) ,
              family = binomial, data = dat1)
summary(model_0)

#-------------------------
# Model with region only
#-------------------------
model_1 <- glmer(cbind(severe_pv_all,sample_size - severe_pv_all) ~   region_1 + (1 | unique_ID) ,
              family = binomial, data = dat1)

summary(model_1)

#--------------------------------------------------
# Model with settings only
#--------------------------------------------------
model_2 <- glmer(cbind(severe_pv_all,sample_size - severe_pv_all) ~   settings_1 + (1 | unique_ID) ,
              family = binomial, data = dat1)

summary(model_2)

#--------------------------------------------------
# Model with region + settings
#--------------------------------------------------

model_3 <- glmer(cbind(severe_pv_all,sample_size - severe_pv_all) ~   region_1 + settings_1 + (1 | unique_ID) ,
              family = binomial, data = dat1)

summary(model_3)

#--------------------------------------------------
# Model with interaction between region+ settings
#--------------------------------------------------

model_4 <- glmer(cbind(severe_pv_all,sample_size - severe_pv_all) ~   region_1*settings_1 + (1 | unique_ID) ,
              family = binomial, data = dat1)
summary(model_4)


# Model comparision using AIC/BIC
anova(model_0,model_1)
anova(model_0,model_2)
anova(model_1,model_3)
anova(model_2,model_3)
anova(model_3,model_4)

anova(model_0,model_1,model_2,model_3,model_4)


#=========================================================================
# Mixed effects logistic regression appraoch: WHO DEFINITION
#=========================================================================

#------------------
# Baseline model
#------------------
model_who_0 <- glmer(cbind(Severe_WHO,sample_size - Severe_WHO) ~    (1 | unique_ID) ,
              family = binomial, data = dat1)
summary(model_who_0)

#-------------------------
# Model with region only
#-------------------------
model_who_1 <- glmer(cbind(Severe_WHO,sample_size - Severe_WHO) ~   region_1 + (1 | unique_ID) ,
              family = binomial, data = dat1)

summary(model_who_1)

#--------------------------------------------------
# Model with settings only
#--------------------------------------------------
model_who_2 <- glmer(cbind(Severe_WHO,sample_size - Severe_WHO) ~   settings_1 + (1 | unique_ID) ,
              family = binomial, data = dat1)

summary(model_who_2)

#--------------------------------------------------
# Model with region + settings
#--------------------------------------------------

model_who_3 <- glmer(cbind(Severe_WHO,sample_size - Severe_WHO) ~   region_1 + settings_1 + (1 | unique_ID) ,
              family = binomial, data = dat1)

summary(model_who_3)

#--------------------------------------------------
# Model with interaction between region+ settings
#--------------------------------------------------

model_who_4 <- glmer(cbind(Severe_WHO,sample_size - Severe_WHO) ~   region_1*settings_1 + (1 | unique_ID) ,
              family = binomial, data = dat1)
summary(model_who_4)


# Model comparision using AIC/BIC
anova(model_who_0,model_who_1)
anova(model_who_0,model_who_2)
anova(model_who_1,model_who_3)
anova(model_who_2,model_who_3)
anova(model_who_3,model_who_4)

# End script
