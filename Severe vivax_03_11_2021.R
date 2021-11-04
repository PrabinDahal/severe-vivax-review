#=================================================================
# Title		: Literature review of severe vivax malaria 
# Data version	: 11-Oct-2021 (from APP)
# Script Date	: 03-Nov-2021
#=================================================================
rm(list=ls())
library(readxl)

# meta-analysis packages
library(meta)
library(dplyr)
library(metasens) # for copas selection

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
dat$cerebral_all 		<- dat$cerebral_WHO+dat$cerebral_other
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
dat1 <- droplevels(dat1 )

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
	dplyr::group_by(region, settings) %>%
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
length(table(dat1$unique_ID))

# Meta-analysis of proportion
(meta.prop <- metaprop(
			data = dat1,
			severe_pv_all,
			sample_size, 
			prediction=TRUE,
			studlab = unique_ID)
		)

# Adjusted for potential small study effect
metabias(meta.prop) 
funnel(meta.prop)
trimfill(meta.prop)
copas(meta.prop)

# Subgroup analysis by region and hopsitalisation settings
update.meta(meta.prop,byvar=region)
update.meta(meta.prop,byvar=settings)

#--------------------------------------------------------
# Adjustment for small study effects separately by region
#--------------------------------------------------------

# Africa only
dat1_afr <- dat1[which(dat1$region=="Africa"),]

(meta.prop.afr <- metaprop(
			data = dat1_afr,
			severe_pv_all,
			sample_size, 
			prediction=TRUE,
			studlab = unique_ID
			)
		)
metabias(meta.prop.afr)
funnel(meta.prop.afr)
trimfill(meta.prop.afr)
copas(meta.prop.afr)
update.meta(meta.prop.afr,byvar=settings)


# Asia only
dat1_asia <- dat1[which(dat1$region=="Asia"),]

(meta.prop.asia<- metaprop(
			data = dat1_asia,
			severe_pv_all,
			sample_size, 
			prediction=TRUE,
			studlab = unique_ID
			)
		)
metabias(meta.prop.asia)
funnel(meta.prop.asia)
trimfill(meta.prop.asia)
copas(meta.prop.asia)
update.meta(meta.prop.asia,byvar=settings)

#-----------------
# South America 
#-----------------
dat1_sa <- dat1[which(dat1$region=="South America"),]

(meta.prop.sa<- metaprop(
			data = dat1_sa,
			severe_pv_all,
			sample_size, 
			prediction=TRUE,
			studlab = unique_ID
			)
		)
metabias(meta.prop.sa)
funnel(meta.prop.sa)
trimfill(meta.prop.sa)
copas(meta.prop.sa)
update.meta(meta.prop.sa,byvar=settings)

#-----------------
# Oceania
#-----------------
table(dat1$region)
dat1_oc <- dat1[which(dat1$region=="Oceania"),]
table(dat1_oc $severe_pv_all)
table(dat1_oc $sample_size)

(meta.prop.oc<- metaprop(
			data = dat1_oc,
			severe_pv_all,
			sample_size, 
			prediction=TRUE,
			studlab = unique_ID
			)
		)
update.meta(meta.prop.oc,byvar=settings)

#-----------------------------------------
# Bias adjustement separately by settings
#-----------------------------------------
#-----------------
# hospitalised
#-----------------
dat1_hosp <- dat1[which(dat1$settings=="hospitalised"),]

(meta.prop.hosp <- metaprop(
			data = dat1_hosp,
			severe_pv_all,
			sample_size, 
			prediction=TRUE,
			studlab = unique_ID
			)
		)
metabias(meta.prop.hosp)
funnel(meta.prop.hosp)
trimfill(meta.prop.hosp)
copas(meta.prop.hosp)

#-----------------
# outpatients
#-----------------
dat1_out <- dat1[which(dat1$settings=="outpatients"),]

(meta.prop.hosp <- metaprop(
			data = dat1_out,
			severe_pv_all,
			sample_size, 
			prediction=TRUE,
			studlab = unique_ID
			)
		)
metabias(meta.prop.hosp)
funnel(meta.prop.hosp)
trimfill(meta.prop.hosp)
copas(meta.prop.hosp)

##################################################################################################
# Meta-analysis II: Estimating proportion of patients with severe vivax based on WHO definition
##################################################################################################
dat2 <- dat[which(!is.na(dat$Severe_WHO) & !is.na(dat$sample_size)),]

dat2 %>% 
	dplyr::summarise(
		n_pub= length(unique(unique_ID)),
		n_severe = sum(Severe_WHO, na.rm=T),
		n_total= sum(sample_size, na.rm=T),
		paste = paste(sum(Severe_WHO, na.rm=T), sum(sample_size, na.rm=T), sep="/")
	)

dat2 %>% 
	dplyr::group_by(region) %>%
	dplyr::summarise(
		n_pub= length(unique(unique_ID)),
		n_severe = sum(Severe_WHO, na.rm=T),
		n_total= sum(sample_size, na.rm=T),
		paste = paste(sum(Severe_WHO, na.rm=T), sum(sample_size, na.rm=T), sep="/")
	)

dat2 %>% 
	dplyr::group_by(settings) %>%
	dplyr::summarise(
		n_pub= length(unique(unique_ID)),
		n_severe = sum(Severe_WHO, na.rm=T),
		n_total= sum(sample_size, na.rm=T),
		paste = paste(sum(Severe_WHO, na.rm=T), sum(sample_size, na.rm=T), sep="/")
	)

dat2 %>% 
	dplyr::group_by(region,settings) %>%
	dplyr::summarise(
		n_pub= length(unique(unique_ID)),
		n_severe = sum(Severe_WHO, na.rm=T),
		n_total= sum(sample_size, na.rm=T),
		paste = paste(sum(Severe_WHO, na.rm=T), sum(sample_size, na.rm=T), sep="/")
	)

(meta.prop.who <- metaprop(
			data = dat2,
			Severe_WHO,
			sample_size, 
			prediction=TRUE,
			studlab = unique_ID
			)
		)

# Adjusted for potential publication biases 
metabias(meta.prop.who)
trimfill(meta.prop.who)
copas(meta.prop.who)


# Subgroup analysis by region & settings
update.meta(meta.prop.who,byvar=region)
update.meta(meta.prop.who,byvar=settings)

#====================
# By regions
#====================

# Africa only
dat2_afr <- dat2[which(dat2$region=="Africa"),]

(meta.prop.who.afr <- metaprop(
			data = dat2_afr,
			Severe_WHO,
			sample_size, 
			prediction=TRUE,
			studlab = unique_ID

			)
		)
metabias(meta.prop.who.afr, k.min=5)
trimfill(meta.prop.who.afr)
copas(meta.prop.who.afr)
update.meta(meta.prop.who.afr,byvar=settings)


# Asia only
dat2_asia <- dat2[which(dat2$region=="Asia"),]

(meta.prop.who.asia<- metaprop(
			data = dat2_asia,
			Severe_WHO,
			sample_size, 
			prediction=TRUE,
			studlab = unique_ID
			)
		)
metabias(meta.prop.who.asia)
trimfill(meta.prop.who.asia)
copas(meta.prop.who.asia)
update.meta(meta.prop.who.asia,byvar=settings)

# South America 
dat2_sa <- dat2[which(dat2$region=="South America"),]

(meta.prop.who.sa<- metaprop(
			data = dat2_sa,
			Severe_WHO,
			sample_size, 
			prediction=TRUE,
			studlab = unique_ID
			)
		)
metabias(meta.prop.who.sa)
trimfill(meta.prop.who.sa)
copas(meta.prop.who.sa)
update.meta(meta.prop.who.sa,byvar=settings)


# South America 
table(dat2$region)
dat2_oc <- dat2[which(dat2$region=="Oceania"),]

(meta.prop.who.oc<- metaprop(
			data = dat2_oc,
			Severe_WHO,
			sample_size, 
			prediction=TRUE,
			studlab = unique_ID
			)
		)
update.meta(meta.prop.who.oc,byvar=settings)


#================
# By settings
#================
dat2_hosp <- dat2[which(dat2$settings=="hospitalised"),]

(meta.prop.who.hosp <- metaprop(
			data = dat2_hosp,
			Severe_WHO,
			sample_size, 
			prediction=TRUE,
			studlab = unique_ID
			)
		)
metabias(meta.prop.who.hosp)
trimfill(meta.prop.who.hosp)
copas(meta.prop.who.hosp)

# outpatients
dat2_out <- dat2[which(dat2$settings=="outpatients"),]

(meta.prop.who.out<- metaprop(
			data = dat2_out,
			Severe_WHO,
			sample_size, 
			prediction=TRUE,
			studlab = unique_ID
			)
		)
metabias(meta.prop.who.out)
trimfill(meta.prop.who.out)
copas(meta.prop.who.out)


##################################################################################################
# Meta-analysis III: Estimating mortality proportion
##################################################################################################
dat4 <- dat[which(!is.na(dat$n_death) & !is.na(dat$sample_size)),]

dat4 %>% 
	dplyr::summarise(
		n_pub= length(unique(unique_ID)),
		n_death = sum(n_death, na.rm=T),
		n_total= sum(sample_size, na.rm=T),
		paste = paste(sum(n_death, na.rm=T), sum(sample_size, na.rm=T), sep="/")
	)

dat4 %>% 
	dplyr::group_by(region) %>%
	dplyr::summarise(
		n_pub= length(unique(unique_ID)),
		n_death = sum(n_death, na.rm=T),
		n_total= sum(sample_size, na.rm=T),
		paste = paste(sum(n_death, na.rm=T), sum(sample_size, na.rm=T), sep="/")
	)

dat4 %>% 
	dplyr::group_by(settings) %>%
	dplyr::summarise(
		n_pub= length(unique(unique_ID)),
		n_death = sum(n_death, na.rm=T),
		n_total= sum(sample_size, na.rm=T),
		paste = paste(sum(n_death, na.rm=T), sum(sample_size, na.rm=T), sep="/")
	)

# overall meta-analysis
(meta.prop.death <- metaprop(
			data = dat4,
			n_death,
			sample_size, 
			prediction=TRUE
			)
		)
# Adjusted for potential publication biases 
metabias(meta.prop.death)
trimfill(meta.prop.death)
copas(meta.prop.death)

#--------------------------------
# Subgroup analysis by region
#--------------------------------
update.meta(meta.prop.death,byvar=region)
update.meta(meta.prop.death,byvar=settings)

#=================
# By Region
#=================

#------------
# Africa
#------------

# There were no events from Africa in 2 studies
# Present simple 95% CI instead

require(binom)
binom.confint(0,2263, method="wilson")

#-----------------
# Asia only
#-----------------
dat4_asia <- dat4[which(dat4$region=="Asia"),]

(meta.prop.mort.asia <- metaprop(
			data = dat4_asia,
			n_death,
			sample_size, 
			prediction=TRUE,
			studlab = unique_ID
			)
		)
metabias(meta.prop.mort.asia)
trimfill(meta.prop.mort.asia)
copas(meta.prop.mort.asia)

#-----------------
# South America 
#-----------------
dat4_sa <- dat4[which(dat4$region=="South America"),]

(meta.prop.mort.sa <- metaprop(
			data = dat4_sa,
			n_death,
			sample_size, 
			prediction=TRUE,
			studlab = unique_ID
			)
		)
metabias(meta.prop.mort.sa)
trimfill(meta.prop.mort.sa)
copas(meta.prop.mort.sa)


#-----------------
# hospitalised
#-----------------
dat4_hosp <- dat4[which(dat4$settings=="hospitalised"),]

(meta.prop.mort.hosp<- metaprop(
			data = dat4_hosp,
			n_death,
			sample_size, 
			prediction=TRUE,
			studlab = unique_ID
			)
		)
metabias(meta.prop.mort.hosp)
trimfill(meta.prop.mort.hosp)
copas(meta.prop.mort.hosp)


#-----------------
# outpatients
#-----------------
dat4_out <- dat4[which(dat4$settings=="outpatients"),]

(meta.prop.mort.out<- metaprop(
			data = dat4_out,
			n_death,
			sample_size, 
			prediction=TRUE,
			studlab = unique_ID
			)
		)
metabias(meta.prop.mort.out)
trimfill(meta.prop.mort.out)
copas(meta.prop.mort.out)

##################################################################################################
# Meta-analysis IV: Estimating proportion of cerebral malaria using WHO
##################################################################################################
dat5 <- dat[which(!is.na(dat$cerebral_WHO) & !is.na(dat$sample_size)),]

dat5 %>% 
	dplyr::summarise(
		n_pub= length(unique(unique_ID)),
		n_cerebral = sum(cerebral_WHO, na.rm=T),
		n_total= sum(sample_size, na.rm=T),
		paste = paste(sum(cerebral_WHO, na.rm=T), sum(sample_size, na.rm=T), sep="/")
	)

dat5 %>% 
	dplyr::group_by(settings) %>%
	dplyr::summarise(
		n_pub= length(unique(unique_ID)),
		n_cerebral = sum(cerebral_WHO, na.rm=T),
		n_total= sum(sample_size, na.rm=T),
		paste = paste(sum(cerebral_WHO, na.rm=T), sum(sample_size, na.rm=T), sep="/")
	)

# meta-analysis
(meta.prop.cerebral.who <- metaprop(
			data = dat5,
			cerebral_WHO,
			sample_size, 
			prediction=TRUE,
			studlab = unique_ID
			)
		)

# assessment and adjustment for publication bias
metabias(meta.prop.cerebral.who)
trimfill(meta.prop.cerebral.who)
copas(meta.prop.cerebral.who)


# sub-group meta-analysis by settings
update.meta(meta.prop.cerebral.who,byvar=settings)

# take further look into hospitalised settings
dat5_hosp <- dat5[which(dat5$settings=="hospitalised"),]

(meta.prop.cerebral.who.hosp <- metaprop(
			data = dat5_hosp,
			cerebral_WHO,
			sample_size, 
			prediction=TRUE,
			studlab = unique_ID
			)
		)
metabias(meta.prop.cerebral.who.hosp)
trimfill(meta.prop.cerebral.who.hosp)
copas(meta.prop.cerebral.who.hosp)


# outpatient settings
dat5_out <- dat5[which(dat5$settings=="outpatients"),]

(meta.prop.cerebral.who.out<- metaprop(
			data = dat5_out,
			cerebral_WHO,
			sample_size, 
			prediction=TRUE,
			studlab = unique_ID
			)
		)
metabias(meta.prop.cerebral.who.out)
trimfill(meta.prop.cerebral.who.out)
copas(meta.prop.cerebral.who.out)


##################################################################################################
# Meta-analysis V: Estimating proportion of cerebral malaria using ALL 
##################################################################################################
dat6 <- dat[which(!is.na(dat$cerebral_all) & !is.na(dat$sample_size)),]

dat6 %>% 
	dplyr::summarise(
		n_pub= length(unique(unique_ID)),
		n_cerebral = sum(cerebral_all, na.rm=T),
		n_total= sum(sample_size, na.rm=T),
		paste = paste(sum(cerebral_all, na.rm=T), sum(sample_size, na.rm=T), sep="/")
	)

dat6 %>% 
	dplyr::group_by(settings) %>%
	dplyr::summarise(
		n_pub= length(unique(unique_ID)),
		n_cerebral = sum(cerebral_all, na.rm=T),
		n_total= sum(sample_size, na.rm=T),
		paste = paste(sum(cerebral_all, na.rm=T), sum(sample_size, na.rm=T), sep="/")
	)

# overall meta-analysis
(meta.prop.cerebral.all  <- metaprop(
			data = dat6,
			cerebral_all,
			sample_size, 
			prediction=TRUE,
			studlab = unique_ID
			)
		)

# testing and adjustment for potential publication bias
metabias(meta.prop.cerebral.all)
trimfill(meta.prop.cerebral.all)
copas(meta.prop.cerebral.all)


# sub-group meta-analysis
update.meta(meta.prop.cerebral.all,byvar=settings)
update.meta(meta.prop.cerebral.all,byvar=region)

# hospitalised settings
dat6_hosp <- dat6[which(dat6$settings=="hospitalised"),]

(meta.prop.cerebral.all.hosp <- metaprop(
			data = dat6_hosp,
			cerebral_all,
			sample_size, 
			prediction=TRUE,
			studlab = unique_ID
			)
		)
metabias(meta.prop.cerebral.all.hosp)
trimfill(meta.prop.cerebral.all.hosp)
copas(meta.prop.cerebral.all.hosp)

# outpatients
dat6_out <- dat6[which(dat6$settings=="outpatients"),]

(meta.prop.cerebral.all.out<- metaprop(
			data = dat6_out,
			cerebral_all,
			sample_size, 
			prediction=TRUE,
			studlab = unique_ID
			)
		)
metabias(meta.prop.cerebral.all.out)
trimfill(meta.prop.cerebral.all.out)
copas(meta.prop.cerebral.all.out)

##################################################################################################
# Meta-analysis VI: Estimating proportion of renal complications using WHO definition
##################################################################################################
dat7 <- dat[which(!is.na(dat$renal_WHO) & !is.na(dat$sample_size)),]

dat7 %>% 
	dplyr::summarise(
		n_pub= length(unique(unique_ID)),
		n_renal = sum(renal_WHO, na.rm=T),
		n_total= sum(sample_size, na.rm=T),
		paste = paste(sum(renal_WHO, na.rm=T), sum(sample_size, na.rm=T), sep="/")
	)

dat7 %>% 
	dplyr::group_by(settings) %>%
	dplyr::summarise(
		n_pub= length(unique(unique_ID)),
		n_renal = sum(renal_WHO, na.rm=T),
		n_total= sum(sample_size, na.rm=T),
		paste = paste(sum(renal_WHO, na.rm=T), sum(sample_size, na.rm=T), sep="/")
	)

(meta.prop.renal.who <- metaprop(
			data = dat7,
			renal_WHO,
			sample_size, 
			prediction=TRUE,
			studlab = unique_ID
			)
		)

metabias(meta.prop.renal.who)
trimfill(meta.prop.renal.who)
copas(meta.prop.renal.who)

# sub-group meta-analysis
update.meta(meta.prop.renal.who,byvar=settings)

# hospitalised settings
dat7_hosp <- dat7[which(dat7$settings=="hospitalised"),]

(meta.prop.renal.who.hosp <- metaprop(
			data = dat7_hosp,
			renal_WHO,
			sample_size, 
			prediction=TRUE,
			studlab = unique_ID
			)
		)
metabias(meta.prop.renal.who.hosp)
trimfill(meta.prop.renal.who.hosp)
copas(meta.prop.renal.who.hosp)

# outpatients settings
dat7_out <- dat7[which(dat7$settings=="outpatients"),]

(meta.prop.renal.who.out <- metaprop(
			data = dat7_out,
			renal_WHO,
			sample_size, 
			prediction=TRUE,
			studlab = unique_ID
			)
		)
metabias(meta.prop.renal.who.out)
trimfill(meta.prop.renal.who.out)
copas(meta.prop.renal.who.out)

##################################################################################################
# Meta-analysis VII: Estimating proportion of renal complications using ALL
##################################################################################################
dat8 <- dat[which(!is.na(dat$renal_all) & !is.na(dat$sample_size)),]

dat8 %>% 
	dplyr::summarise(
		n_pub= length(unique(unique_ID)),
		n_renal = sum(renal_all, na.rm=T),
		n_total= sum(sample_size, na.rm=T),
		paste = paste(sum(renal_all, na.rm=T), sum(sample_size, na.rm=T), sep="/")
	)

dat8 %>% 
	dplyr::group_by(settings) %>%
	dplyr::summarise(
		n_pub= length(unique(unique_ID)),
		n_renal = sum(renal_all, na.rm=T),
		n_total= sum(sample_size, na.rm=T),
		paste = paste(sum(renal_all, na.rm=T), sum(sample_size, na.rm=T), sep="/")
	)

(meta.prop.renal.all <- metaprop(
			data = dat8,
			renal_all,
			sample_size, 
			prediction=TRUE,
			studlab = unique_ID

			)
		)

# testing and adjustment for potential publication biases
metabias(meta.prop.renal.all)
trimfill(meta.prop.renal.all)
copas(meta.prop.renal.all)

# sub-group meta-analysis
update.meta(meta.prop.renal.all,byvar=settings)

# hospitalised settings
dat8_hosp <- dat8[which(dat8$settings=="hospitalised"),]

(meta.prop.renal.all.hosp <- metaprop(
			data = dat8_hosp,
			renal_all,
			sample_size, 
			prediction=TRUE,
			studlab = unique_ID
			)
		)
metabias(meta.prop.renal.all.hosp)
trimfill(meta.prop.renal.all.hosp)
copas(meta.prop.renal.all.hosp)

# outpatients settings
dat8_out <- dat8[which(dat8$settings=="outpatients"),]

(meta.prop.renal.all.out<- metaprop(
			data = dat8_out,
			renal_all,
			sample_size, 
			prediction=TRUE,
			studlab = unique_ID
			)
		)
metabias(meta.prop.renal.all.out)
trimfill(meta.prop.renal.all.out)
copas(meta.prop.renal.all.out)


##################################################################################################
# Meta-analysis IX: Estimating proportion of respitory complications using WHO definition
##################################################################################################
dat9 <- dat[which(!is.na(dat$respiratory_WHO) & !is.na(dat$sample_size)),]

dat9 %>% 
	dplyr::summarise(
		n_pub= length(unique(unique_ID)),
		n_resp = sum(respiratory_WHO, na.rm=T),
		n_total= sum(sample_size, na.rm=T),
		paste = paste(sum(respiratory_WHO, na.rm=T), sum(sample_size, na.rm=T), sep="/")
	)

dat9 %>% 
	dplyr::group_by(settings) %>%
	dplyr::summarise(
		n_pub= length(unique(unique_ID)),
		n_resp = sum(respiratory_WHO, na.rm=T),
		n_total= sum(sample_size, na.rm=T),
		paste = paste(sum(respiratory_WHO, na.rm=T), sum(sample_size, na.rm=T), sep="/")
	)

(meta.prop.resp.who <- metaprop(
			data = dat9,
			respiratory_WHO,
			sample_size, 
			prediction=TRUE,
			studlab = unique_ID
			)
		)

metabias(meta.prop.resp.who)
trimfill(meta.prop.resp.who)
copas(meta.prop.resp.who)

# sub-group meta-analysis
update.meta(meta.prop.resp.who,byvar=settings)

# hospitalised settings
dat9_hosp <- dat9[which(dat9$settings=="hospitalised"),]

(meta.prop.resp.who.hosp <- metaprop(
			data = dat9_hosp,
			respiratory_WHO,
			sample_size, 
			prediction=TRUE,
			studlab = unique_ID
			)
		)
metabias(meta.prop.resp.who.hosp)
trimfill(meta.prop.resp.who.hosp)
copas(meta.prop.resp.who.hosp)

# otupatients settings
dat9_out  <- dat9[which(dat9$settings=="outpatients"),]

(meta.prop.resp.who.out <- metaprop(
			data = dat9_out,
			respiratory_WHO,
			sample_size, 
			prediction=TRUE,
			studlab = unique_ID
			)
		)
metabias(meta.prop.resp.who.out)
trimfill(meta.prop.resp.who.out)
copas(meta.prop.resp.who.out)

##################################################################################################
# Meta-analysis X: Estimating proportion of respitory complications using ALL definition
##################################################################################################
dat10 <- dat[which(!is.na(dat$respiratory_all) & !is.na(dat$sample_size)),]

dat10 %>% 
	dplyr::summarise(
		n_pub= length(unique(unique_ID)),
		n_resp = sum(respiratory_all, na.rm=T),
		n_total= sum(sample_size, na.rm=T),
		paste = paste(sum(respiratory_all, na.rm=T), sum(sample_size, na.rm=T), sep="/")
	)

dat10 %>% 
	dplyr::group_by(settings) %>%
	dplyr::summarise(
		n_pub= length(unique(unique_ID)),
		n_resp = sum(respiratory_all, na.rm=T),
		n_total= sum(sample_size, na.rm=T),
		paste = paste(sum(respiratory_all, na.rm=T), sum(sample_size, na.rm=T), sep="/")
	)

# overall meta-analysis
(meta.prop.resp.all<- metaprop(
			data = dat10,
			respiratory_all,
			sample_size, 
			prediction=TRUE,
			studlab = unique_ID
			)
		)

# assessment and adjustment for potential publication bias
metabias(meta.prop.resp.all)
trimfill(meta.prop.resp.all)
copas(meta.prop.resp.all)

# sub-group meta-analysis
update.meta(meta.prop.resp.all,byvar=settings)

# hospitalised settings
dat10_hosp <- dat10[which(dat10$settings=="hospitalised"),]

(meta.prop.resp.all.hosp <- metaprop(
			data = dat10_hosp,
			respiratory_all,
			sample_size, 
			prediction=TRUE,
			studlab = unique_ID
			)
		)
metabias(meta.prop.resp.all.hosp)
trimfill(meta.prop.resp.all.hosp)
copas(meta.prop.resp.all.hosp)


# outpatients settings
dat10_out  <- dat10[which(dat10$settings=="outpatients"),]

(meta.prop.resp.all.out <- metaprop(
			data = dat10_out,
			respiratory_all,
			sample_size, 
			prediction=TRUE,
			studlab = unique_ID
			)
		)
metabias(meta.prop.resp.all.out)
trimfill(meta.prop.resp.all.out)
copas(meta.prop.resp.all.out)

#########################################################################################
# Meta-analysis XI: Estimating proportion of patients with severe vivax in pregnancy
#########################################################################################

# Severe vivax
preg %>% 
	dplyr::summarise(
		n_pub= length(unique(unique_ID)),
		n_sev_pv = sum(severe_pv_all, na.rm=T),
		n_total= sum(sample_size, na.rm=T),
		paste = paste(sum(severe_pv_all, na.rm=T), sum(sample_size, na.rm=T), sep="/")
	)

preg %>% 
	dplyr::summarise(
		n_pub= length(unique(unique_ID)),
		n_sev_pv = sum(Severe_WHO, na.rm=T),
		n_total= sum(sample_size, na.rm=T),
		paste = paste(sum(Severe_WHO, na.rm=T), sum(sample_size, na.rm=T), sep="/")
	)

# Mortality
preg %>% 
	dplyr::summarise(
		n_pub= length(unique(unique_ID)),
		n_death = sum(n_death, na.rm=T),
		n_total= sum(sample_size, na.rm=T),
		paste = paste(sum(n_death, na.rm=T), sum(sample_size, na.rm=T), sep="/")
	)

# Respiratory
preg %>% 
	dplyr::summarise(
		n_pub= length(unique(unique_ID)),
		n_death = sum(respiratory_all, na.rm=T),
		n_total= sum(sample_size, na.rm=T),
		paste = paste(sum(respiratory_all, na.rm=T), sum(sample_size, na.rm=T), sep="/")
	)

preg %>% 
	dplyr::summarise(
		n_pub= length(unique(unique_ID)),
		n_death = sum(respiratory_WHO, na.rm=T),
		n_total= sum(sample_size, na.rm=T),
		paste = paste(sum(respiratory_WHO, na.rm=T), sum(sample_size, na.rm=T), sep="/")
	)

# Cerebral
preg %>% 
	dplyr::summarise(
		n_pub= length(unique(unique_ID)),
		n_cerebral = sum(cerebral_all, na.rm=T),
		n_total= sum(sample_size, na.rm=T),
		paste = paste(sum(cerebral_all, na.rm=T), sum(sample_size, na.rm=T), sep="/")
	)

preg %>% 
	dplyr::summarise(
		n_pub= length(unique(unique_ID)),
		n_cerebral = sum(cerebral_WHO, na.rm=T),
		n_total= sum(sample_size, na.rm=T),
		paste = paste(sum(cerebral_WHO, na.rm=T), sum(sample_size, na.rm=T), sep="/")
	)

# Renal
preg %>% 
	dplyr::summarise(
		n_pub= length(unique(unique_ID)),
		n_renal = sum(renal_all, na.rm=T),
		n_total= sum(sample_size, na.rm=T),
		paste = paste(sum(renal_all, na.rm=T), sum(sample_size, na.rm=T), sep="/")
	)

preg %>% 
	dplyr::summarise(
		n_pub= length(unique(unique_ID)),
		n_renal = sum(renal_WHO, na.rm=T),
		n_total= sum(sample_size, na.rm=T),
		paste = paste(sum(renal_WHO, na.rm=T), sum(sample_size, na.rm=T), sep="/")
	)

# proportion of severe vivax malaria: all
(meta.prop.preg  <- metaprop(
			data = preg ,
			severe_pv_all,
			sample_size, 
			prediction=TRUE,
			studlab = unique_ID
			)
		)

## proportion of severe vivax malaria: WHO

(meta.prop.preg.who  <- metaprop(
			data = preg ,
			Severe_WHO,
			sample_size, 
			prediction=TRUE,
			studlab = unique_ID
			)
		)

## Mortality
(meta.prop.preg.mort <- metaprop(
			data = preg ,
			n_death,
			sample_size, 
			prediction=TRUE,
			studlab = unique_ID
			)
		)

# End Script(Not Run)