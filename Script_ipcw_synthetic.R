
# File name: Script_ipcw_synthetic.R
# Author: Amani Al Tawil
# Date: March 7, 2024
# Data set: df_ipcw_pp (synthetic data in person-time format) 
# Method: IPCW analyses to adjust for Loss to Follow-up (LTFU)/Administrative Censoring(AC) 
#         &/or Treatment Switching (XO) 
# References: Robins 2000 - Cole and Herna'n 2008 - Sullivan 2020 - Murray 2021
#
#
#
# Description:
# 
# This R script performs Inverse Probability of Censoring Weighting (IPCW) analyses
# on the generated synthetic data "df_ipcw_pp" 
# to estimate treatment effects while eliminating LTFU/AC and/or XO. 
# This analysis follows the methodology described by Robins & Finkelstein 2000. 
# Refer to Electronic Supplementary Material (ESM) 1 for the data dictionary.
#
#
#
# Notes: 
# 
# To access the code that is used to generate the synthetic data "df_ipcw_pp", 
# refer to the R script 'synthetic_data.R'.
# This R script should not be used to 
# reproduce the published results or draw clinical conclusions. 
# This script should only serve as a tool for understanding our analysis process.
# For details about the variables included in the weighting models, refer to ESM1, 
# ESM 3 and Tables 4 and 6 in the manuscript.
# Since we only provide this code applied to the synthetic data,
# we do not include details regarding our model selection process (irrelevant).
# Running these analyses here might yield extreme weights, or unexpected results
# This is because the variables included in the weighting models are selected 
# based on our original analyses
# To manage excessive computational time and for demonstration purposes,
# we have set the number of bootstrap samples (R= 20) when computing CI.
# For your analyses, consider increasing R to 1000.
#
#
#
# Analysis 1: LTFU/AC Analysis:
#
# In lines 184 to 718
# We perform six analyses to estimate the ITT effect of the two treatment strategies as follows:
# “Assign to crizotinib” (Z = 0, $\bar{C}$ = 0): This strategy assigns crizotinib at baseline and enforces no LTFU vs
# “Assign to brigatinib” (Z = 1, $\bar{C}$ = 0): This strategy assigns brigatinib at baseline and enforces no LTFU
# This effect is investigated through various models (specification 1 to 6) for the construction of weights for LTFU/AC;
# Death censored by LTFU/AC only
#
#
#
# Analysis 2: Switching Analysis:
#
# In lines 723 to 1470 
# We perform eight analyses to estimate the causal effect of
# “Always treat with crizotinib” ($\bar{A}$ = 0): This strategy administers crizotinib at each time point vs
# “Always treat with brigatinib” ($\bar{A}$ = 1): This strategy administers brigatinib at each time point 
# i.e. the effect of hypothetical treatment strategies on overall survival that do not allow for treatment switching
# This effect is investigated through various models (specification 1 to 8) for the construction of weights for switching;
# LTFU assumed at random; Deaths censored by a minimum of treatment switching and LTFU/AC
# In lines 1455 to 1470 We apply the analyses from spec 7 with truncated weights
#
#
#
# Analysis 3: LTFU/AC and Switching Analysis:
#
# In lines 1472 to 2288 
# We perform two analyses that adjusts for both XO and LTFU to estimate the causal effect of
# “Always treat with crizotinib” ($\bar{A}$ = 0, $\bar{C}$ = 0): This strategy administers crizotinib at each time point and enforces no LTFU
# “Always treat with brigatinib” ($\bar{A}$ = 1, $\bar{C}$ = 0): This strategy administers brigatinib at each time point and enforces no LTFU
# This effect is investigated through models (specification 4 and 7) for the construction of weights for switching and LTFU/AC;
# multiplication of XO and LTFU weights is applied in the outcome model
# Deaths censored by a minimum of treatment switching and LTFU/AC
#
#
#
# Usage:
#
# 1. Set the working directory to the location of your data set and this script.
# 2. Ensure that the required libraries are installed and loaded.
# 3. Customize the parameters and variables as needed for your specific analysis.
# 4. Run the script to perform the IPCW analysis.
#
#
#
# Set Working directory---------------------------------------------------------

# Replace "YOUR_DIRECTORY_PATH" with the actual directory path.
setwd("YOUR_DIRECTORY_PATH") 
#setwd("D:\\LMU--2020\\PhD\\BMC\\final_submission\\switching_github")
getwd()
dir()
#
#
#
# Load data sets----------------------------------------------------------------

options("install.lock"=FALSE)
load("df_ipcw_pp.Rda")


# Prevent scientific notation
options(scipen=999)


# In our actual analyses we set K=48 months
# Here we set k to 36 since our synthetic data contains only 3 observation with 48 months
# Set k (for time) to 36 months
k<- 36
#
#
#
# Install required Packages-----------------------------------------------------

#This will automatically install them into your default library
if(!require("DescTools"))install.packages("DescTools")
library(DescTools)
if(!require("dplyr"))install.packages("dplyr")
library(dplyr)
if(!require("survival"))install.packages("survival")
library(survival)
if(!require("survminer"))install.packages("survminer")
library(survminer)
if(!require("ggplot2"))install.packages("ggplot2")
library(ggplot2)
if(!require("boot"))install.packages("boot")
library(boot)
if(!require("Hmisc"))install.packages("Hmisc")
library(Hmisc)
if(!require("lmtest"))install.packages("lmtest")
library(lmtest)
if(!require("sandwich"))install.packages("sandwich")
library(sandwich)
if(!require("miceadds"))install.packages("miceadds")
library(miceadds)
if(!require("MASS"))install.packages("MASS")
library(MASS)
if(!require("table1"))install.packages("table1")
library(table1)
if(!require("openxlsx"))install.packages("openxlsx")
library(openxlsx)
if(!require("gee"))install.packages("gee")
library(gee)
if(!require("geepack"))install.packages("geepack")
library(geepack)
if(!require("stats"))install.packages("stats")
library(stats)
if(!require("data.table"))install.packages("data.table")
library(data.table)
if(!require("forestplot"))install.packages("forestplot")
library(forestplot)
if(!require("rms"))install.packages("rms")
library(rms)
#
#
#
# Helper function which generates bootstrap data--------------------------------

create_bootstrap_data <- function(data){
  # Selecting which IDs are included in the bootstrap sample
  n_ids <- length(unique(data$id))
  ids_dt <- as.data.table(sample(1:n_ids, n_ids, replace = TRUE))
  ids_dt[, 'bid' := 1:n_ids]
  colnames(ids_dt) <- c("id", "bid")
  
  # Merging the new ID data set to the original data set
  resample_data <- as.data.table(data)
  setkey(resample_data, "id")
  resample_data <- resample_data[J(ids_dt), allow.cartesian = TRUE]
  resample_data[, 'id' := resample_data$bid]
  resample_data[, 'bid' := NULL]
  
  return(as.data.frame(resample_data))
}
#
#
#
#-------------------------------------------------------------------------------




# Analysis 1: Lost to follow-up/Administrative censoring adjusted analysis------------------
# Analysis process for Table 6 in manuscript
#
#
#
# Commenting the code:
# In lines 223 to 450 entitled "Spec 1 to Spec 6"
# We provide the full code for the 6 analyses with explanatory comments 
# uncomment and execute only one specification at a time
# For demonstration we uncomment spec 1 and comment specs 2 to 6
# All subsequent steps (weight calculation and the different outcome models) should be run with any spec
# In lines 682 to 718 and for demonstration we provide the code for the forest plot with the results from the real data
#
#
#
# Load data
df_censoring<- df_ipcw_pp

analyses.boot<- function(data, indices){
  
  ## Step 1: Generated bootstrap data, dd
  
  if (identical(indices, 1:length(indices))){
    dd <- data
  } else {
    dd <- create_bootstrap_data(data)
  }
  
  dd$wt1 <- NULL
  dd$stabwt1<- NULL
  # Run the below lines for spec 3 and 5
  knots0 = Hmisc::rcspline.eval(dd$time, knots.only = T, pc=T)
  knots1 = Hmisc::rcspline.eval(dd$AGE, knots.only = T, pc=T)
  knots2 = Hmisc::rcspline.eval(df_censoring$tltd_wins, knots.only = T, pc=T)
  
  ## Step 2: Fit the weight models based on dd and calculate the weights
  ### Important: uncomment and execute only one specification at a time
  ### For demonstration we uncomment spec 1 and comment specs 2 to 6
  
  ### Spec 1--------------------------------------------------------------------

  # Numerator control for LTFU
  nmodel_lofu_0<- glm(censoringtd_new ~ time + time2,
                       data=dd[dd$plannedtrt==0,], family =binomial())

  # Numerator exp for LTFU
   nmodel_lofu_1<- glm(censoringtd_new ~ time + time2,
                      data=dd[dd$plannedtrt==1,], family =binomial())

  # Predict
   dd$pnumerator1_0<- predict(nmodel_lofu_0, newdata=dd, type='response')
   dd$pnumerator1_1<- predict(nmodel_lofu_1, newdata=dd, type='response')

  # Denominator control for LTFU
   dmodel_lofu_0<- glm(
     censoringtd_new ~ time + SMKHISGR1 + STRATRAN + ecogtdGR2
     + AGE + MICNSFL + tltd_wins + progtd + trttd + LIASE + DIAGINITGR3
    + time2
     + AGE2,
     data=dd[dd$plannedtrt==0,],
     family =binomial())

  # Denominator exp for LTFU
   dmodel_lofu_1<- glm(
     censoringtd_new ~ time + AGE + SEX + LIASE + STRATRAN + ECOGGR1 + ecogtdGR2
     + tltd_wins + icprogtd + progtd + RACEGR1 + SMKHISGR1 + MICNSFL
     + DIAGINITGR3 + PRRADYN + time2 + AGE2,
     data=dd[dd$plannedtrt==1,],
     family =binomial())

  # Predict
  dd$pdenominator1_0<- predict(dmodel_lofu_0, newdata=dd, type='response')
  dd$pdenominator1_1<- predict(dmodel_lofu_1, newdata=dd, type='response')
  #
  # 
  ### Spec 2----------------------------------------------------------------------

  # # Numerator control for LTFU
  # nmodel_lofu_0<- glm(censoringtd_new ~ time + time2,
  #                     data=dd[dd$plannedtrt==0,], family =binomial())
  # 
  # # Numerator exp for LTFU
  # nmodel_lofu_1<- glm(censoringtd_new ~ time + time2,
  #                     data=dd[dd$plannedtrt==1,], family =binomial())
  # 
  # # Predict
  # dd$pnumerator1_0<- predict(nmodel_lofu_0, newdata=dd, type='response')
  # dd$pnumerator1_1<- predict(nmodel_lofu_1, newdata=dd, type='response')
  # 
  # # Denominator control for LTFU
  # dmodel_lofu_0<- glm(
  #   censoringtd_new ~ time  + STRATRAN
  #   + AGE + DIAGINITGR3 + tltd_wins
  #   + time2
  #   + AGE2
  #   + progtd,
  #   data=dd[dd$plannedtrt==0,],
  #   family =binomial())
  # 
  # # Denominator exp for LTFU
  # dmodel_lofu_1<- glm(
  #   censoringtd_new ~ time + SEX + ECOGGR1 + ecogtdGR2 + icprogtd + progtd + SMKHISGR1
  #   + time2
  #   + DIAGINITGR3
  #   + AGE + AGE2,
  #   data=dd[dd$plannedtrt==1,],
  #   family =binomial())
  # 
  # # Predict
  # dd$pdenominator1_0<- predict(dmodel_lofu_0, newdata=dd, type='response')
  # dd$pdenominator1_1<- predict(dmodel_lofu_1, newdata=dd, type='response')
  # 
  # 
  ### Spec 3----------------------------------------------------------------------
# 
#   # Numerator control for LTFU
#   nmodel_lofu_0<- glm(censoringtd_new ~ rcs(time,knots0),
#                       data=dd[dd$plannedtrt==0,], family =binomial())
# 
#   # Numerator exp for LTFU
#   nmodel_lofu_1<- glm(censoringtd_new ~ rcs(time,knots0),
#                       data=dd[dd$plannedtrt==1,], family =binomial())
# 
#   # Predict
#   dd$pnumerator1_0<- predict(nmodel_lofu_0, newdata=dd, type='response')
#   dd$pnumerator1_1<- predict(nmodel_lofu_1, newdata=dd, type='response')
# 
#   # Denominator control for LTFU
#   dmodel_lofu_0<- glm(
#     censoringtd_new ~ rcs(time,knots0)
#     + SMKHISGR1 + STRATRAN + ecogtdGR2
#     + rcs(AGE,knots1)
#     + MICNSFL
#     + rcs(tltd_wins,knots2)
#     + progtd + trttd + LIASE + DIAGINITGR3,
#     data=dd[dd$plannedtrt==0,],
#     family =binomial())
# 
#   # Denominator exp for LTFU
#   dmodel_lofu_1<- glm(
#     censoringtd_new ~ rcs(time,knots0)
#     + rcs(AGE,knots1)
#     + SEX + LIASE + STRATRAN + ECOGGR1 + ecogtdGR2
#     + rcs(tltd_wins,knots2)
#     + icprogtd + progtd + RACEGR1 + SMKHISGR1 + MICNSFL
#     + DIAGINITGR3 + PRRADYN,
#     data=dd[dd$plannedtrt==1,],
#     family =binomial())
# 
#   # Predict
#   dd$pdenominator1_0<- predict(dmodel_lofu_0, newdata=dd, type='response')
#   dd$pdenominator1_1<- predict(dmodel_lofu_1, newdata=dd, type='response')
  # 
  # 
  ### Spec 4----------------------------------------------------------------------

  # # Numerator control for LTFU
  # nmodel_lofu_0<- glm(censoringtd_new ~ time + time2,
  #                     data=dd[dd$plannedtrt==0,], family =binomial())
  # 
  # # Numerator exp for LTFU
  # nmodel_lofu_1<- glm(censoringtd_new ~ time + time2,
  #                     data=dd[dd$plannedtrt==1,], family =binomial())
  # 
  # # Predict
  # dd$pnumerator1_0<- predict(nmodel_lofu_0, newdata=dd, type='response')
  # dd$pnumerator1_1<- predict(nmodel_lofu_1, newdata=dd, type='response')
  # 
  # # Denominator control for LTFU
  # dmodel_lofu_0<- glm(
  #   censoringtd_new ~ time + SMKHISGR1 + STRATRAN + ecogtdGR1
  #   + AGE + MICNSFL + tltd_wins + progtd + trttd + LIASE + DIAGINITGR3
  #   + time2
  #   + AGE2,
  #   data=dd[dd$plannedtrt==0,],
  #   family =binomial())
  # 
  # # Denominator exp for LTFU
  # dmodel_lofu_1<- glm(
  #   censoringtd_new ~ time + AGE + SEX + LIASE + STRATRAN + ECOGGR1 + ecogtdGR1
  #   + tltd_wins + icprogtd + progtd + RACEGR1 + SMKHISGR1 + MICNSFL
  #   + DIAGINITGR3 + PRRADYN + time2 + AGE2,
  #   data=dd[dd$plannedtrt==1,],
  #   family =binomial())
  # 
  # # Predict
  # dd$pdenominator1_0<- predict(dmodel_lofu_0, newdata=dd, type='response')
  # dd$pdenominator1_1<- predict(dmodel_lofu_1, newdata=dd, type='response')

  # 
  ### Spec 5----------------------------------------------------------------------
  # 
  # # Numerator control for LTFU
  # nmodel_lofu_0<- glm(censoringtd_new ~ rcs(time,knots0),
  #                     data=dd[dd$plannedtrt==0,], family =binomial())
  # 
  # # Numerator exp for LTFU
  # nmodel_lofu_1<- glm(censoringtd_new ~ rcs(time,knots0),
  #                     data=dd[dd$plannedtrt==1,], family =binomial())
  # 
  # # Predict
  # dd$pnumerator1_0<- predict(nmodel_lofu_0, newdata=dd, type='response')
  # dd$pnumerator1_1<- predict(nmodel_lofu_1, newdata=dd, type='response')
  # 
  # # Denominator control for LTFU
  # dmodel_lofu_0<- glm(
  #   censoringtd_new ~ rcs(time,knots0)
  #   + STRATRAN
  #   + rcs(AGE,knots1)
  #   + DIAGINITGR3
  #   + rcs(tltd_wins,knots2)
  #   + progtd,
  #   data=dd[dd$plannedtrt==0,],
  #   family =binomial())
  # 
  # # Denominator exp for LTFU
  # dmodel_lofu_1<- glm(
  #   censoringtd_new ~ rcs(time,knots0)
  #   + SEX + ECOGGR1 + ecogtdGR2 + icprogtd + progtd + SMKHISGR1
  #   + DIAGINITGR3
  #   + rcs(AGE,knots1),
  #   data=dd[dd$plannedtrt==1,],
  #   family =binomial())
  # 
  # # Predict
  # dd$pdenominator1_0<- predict(dmodel_lofu_0, newdata=dd, type='response')
  # dd$pdenominator1_1<- predict(dmodel_lofu_1, newdata=dd, type='response')

  # 
  ### Spec 6----------------------------------------------------------------------
  # 
  # # Numerator control for LTFU
  # nmodel_lofu_0<- glm(censoringtd_new ~ time + time2,
  #                     data=dd[dd$plannedtrt==0,], family =binomial())
  # 
  # # Numerator exp for LTFU
  # nmodel_lofu_1<- glm(censoringtd_new ~ time + time2,
  #                     data=dd[dd$plannedtrt==1,], family =binomial())
  # 
  # # Predict
  # dd$pnumerator1_0<- predict(nmodel_lofu_0, newdata=dd, type='response')
  # dd$pnumerator1_1<- predict(nmodel_lofu_1, newdata=dd, type='response')
  # 
  # # Denominator control for LTFU
  # dmodel_lofu_0<- glm(
  #   censoringtd_new ~ time  + STRATRAN
  #   + AGE + DIAGINITGR3 + tltd_wins
  #   + time2
  #   + AGE2
  #   + progtd,
  #   data=dd[dd$plannedtrt==0,],
  #   family =binomial())
  # 
  # # Denominator exp for LTFU
  # dmodel_lofu_1<- glm(
  #   censoringtd_new ~ time + AGE + ecogtdGR2 + icprogtd + progtd + SMKHISGR1
  #   + time2
  #   + AGE2,
  #   data=dd[dd$plannedtrt==1,],
  #   family =binomial())
  # 
  # # Predict
  # dd$pdenominator1_0<- predict(dmodel_lofu_0, newdata=dd, type='response')
  # dd$pdenominator1_1<- predict(dmodel_lofu_1, newdata=dd, type='response')
  # 
  # 
  #
  ### Weight Calculation----------------------------------------------------------
  
  # Order the data by "id" and "time"
  dd<- dd[order(dd$id,dd$time),]
  
  # Calculate num1 and den1 
  dd<- dd%>%dplyr::mutate(
    num1= ifelse(plannedtrt==1,
                 (censoringtd_new*pnumerator1_1 + (1-censoringtd_new)*(1-pnumerator1_1)),
                 (censoringtd_new*pnumerator1_0 + (1-censoringtd_new)* (1-pnumerator1_0))),
    den1=ifelse(plannedtrt==1,
                (censoringtd_new*pdenominator1_1 + (1-censoringtd_new)* (1-pdenominator1_1)),
                (censoringtd_new*pdenominator1_0 + (1-censoringtd_new)* (1-pdenominator1_0))),
    num1= ifelse(time==0,1,num1),
    den1= ifelse(time==0,1,den1)
  ) %>%
    dplyr::group_by(id) %>%
    dplyr::mutate(
      n1= cumprod(num1),
      d1= cumprod(den1)
    )%>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      stabwt1= n1/d1, # stabilized weights
      wt1=1/d1        # unstabilized weights 
    )
  
  dd_wts= dd[dd$censoringtd_new==0,]
  
  # Distribution of wts
  mean_unstab_ltfu<- mean(dd_wts$wt1)
  sd_unstab_ltfu<- sd(dd_wts$wt1)
  min_unstab_ltfu<- min(dd_wts$wt1)
  max_unstab_ltfu<- max(dd_wts$wt1)
  
  
  # Distribution of Stabwts
  mean_stab_ltfu<- mean(dd_wts$stabwt1)
  sd_stab_ltfu<- sd(dd_wts$stabwt1)
  min_stab_ltfu<- min(dd_wts$stabwt1)
  max_stab_ltfu<- max(dd_wts$stabwt1)
  
  
  
  ## Step 3: Remove observations at censoring and fit the outcomes model--------
  ## replace weights= stabwt1 with weights= wt1 for analyses with unstabilized weights
  
  dd_outcome<- dd[dd$censoringtd_new==0,]

  # Cox proportional hazard model

  ltfu.cox <- coxph(Surv(time,time+1,deathtd) ~ plannedtrt + cluster(id), robust = T,
                   data = dd_outcome, weights= stabwt1, ties = ("breslow")) # weights= wt1
  
  ### HR
  HR<- exp(coef(ltfu.cox))[1]
  
  
  # Pooled logistic regression (dHR, cHR, RR)
  ltfu.logit<- glm(deathtd ~ time + time2 + factor(plannedtrt),
    data= dd_outcome, weights= stabwt1, family =quasibinomial()) # weights= wt1

  ### dHR
  dHR<- exp(coef(ltfu.logit))[4]
  ### cHR and RR
  dd_outcome$prob<- 1-predict(ltfu.logit, newdata=dd_outcome, type="response") # prob: survival probability
  dd_1<- dd_outcome%>%
    dplyr::arrange(id,time)%>%
    dplyr::group_by(id)%>%
    dplyr::mutate(
      surv=cumprod(prob))
  dd_1<- dd_1[,c('surv','plannedtrt','time')]
  results<- dd_1%>%
    dplyr::group_by(time, plannedtrt)%>%
    dplyr::summarize(mean_survival=mean(surv))
  results<- results%>%
    dplyr::ungroup()%>%
    dplyr:: mutate(
      time=time+1
    )
  results<- dplyr::bind_rows(c(time=0,plannedtrt=0,mean_survival=1),c(time=0,plannedtrt=1,mean_survival=1),results )
  results$plannedtrtf<- factor(results$plannedtrt,labels=c("Crizotinib","Brigatinib"))
  wide_results<- reshape2::dcast(results, time ~ plannedtrtf,value.var = 'mean_survival')
  wide_results<- wide_results%>%
    mutate(
      risk0=(1-Crizotinib),
      risk1=(1-Brigatinib),
      RD= (1-Brigatinib)-(1-Crizotinib),
      cHR= log(Brigatinib)/log(Crizotinib),
      RR=(1-Brigatinib)/(1-Crizotinib))


  # KM Estimator (Note: run K= 36)
  kmfit<- survfit(Surv(time, time+1,deathtd) ~ plannedtrt, weights = stabwt1, data = dd_outcome) # weights= wt1

  ### cHR
  S0<- summary(kmfit, times=k, extend=T)$surv[1]
  S1<- summary(kmfit, times=k, extend=T)$surv[2]
  cHR_km<- log(S1)/log(S0)
  # 
  ### RR
  R0_km<- 1-summary(kmfit, times = k)$surv[1]
  R1_km<- 1-summary(kmfit, times = k)$surv[2]
  RD_km<- R1_km-R0_km
  RR_km<- R1_km/R0_km
  
  
  
  # Return results--------------------------------------------------------------
  
  return(c(HR,
           dHR,
           wide_results$risk0[which(wide_results$time==k-1)],
           wide_results$risk1[which(wide_results$time==k-1)],
           wide_results$RD[which(wide_results$time==k-1)],
           wide_results$cHR[which(wide_results$time==k-1)],
           wide_results$RR[which(wide_results$time==k-1)],
           cHR_km,
           RR_km,
           mean_unstab_ltfu,
           sd_unstab_ltfu,
           min_unstab_ltfu,
           max_unstab_ltfu,
           mean_stab_ltfu,
           sd_stab_ltfu,
           min_stab_ltfu,
           max_stab_ltfu))
} 

# Set replicate failures to NA
boot.robust<- function(data,indices){
  tryCatch({analyses.boot(data,indices)}, error= function(e){NA})
}

set.seed(7268)

final.results<- boot(data= df_censoring, statistic= boot.robust, R=20) # 1000 
final.results$t0

# CoxPH
HR<- round(final.results$t0[[1]],2)
HR
conf_interval_HR <- boot.ci(final.results, conf = 0.95, type = "perc", index = 1)
conf_interval_HR


# Pooled logistic regression
dHR<- round(final.results$t0[[2]],2)
dHR
conf_interval_dHR <- boot.ci(final.results, conf = 0.95, type = "perc", index = 2)
conf_interval_dHR
#
cHR<- round(final.results$t0[[6]],2)
cHR
conf_interval_cHR <- boot.ci(final.results, conf = 0.95, type = "perc", index = 6)
conf_interval_cHR
#
RR<- round(final.results$t0[[7]],2)
RR
conf_interval_RR <- boot.ci(final.results, conf = 0.95, type = "perc", index = 7)
conf_interval_RR


# KM Estimator
cHR_km<-round(final.results$t0[[8]],2)
cHR_km
conf_interval_cHR_km<- boot.ci(final.results, conf = 0.95, type = "perc", index = 8)
conf_interval_cHR_km
#
RR_km<-round(final.results$t0[[9]],2)
RR_km
conf_interval_RR_km<- boot.ci(final.results, conf = 0.95, type = "perc", index = 9)
conf_interval_RR_km
#
#
#
# Distribution of LTFU unstab wts
mean_unstab_ltfu<- round(as.numeric(final.results$t0[[10]]),2)
sd_unstab_ltfu<- round(as.numeric(final.results$t0[[11]]),2)
mean_sd_unstab_ltfu<- paste0(mean_unstab_ltfu,"(",sd_unstab_ltfu,")")
min_unstab_ltfu<- round(final.results$t0[[12]],2)
max_unstab_ltfu<- round(final.results$t0[[13]],2)
min_max_unstab_ltfu<- paste0(min_unstab_ltfu,"(",max_unstab_ltfu,")")

# Distribution of LTFU stab wts
mean_stab_ltfu<- round(as.numeric(final.results$t0[[14]]),2)
sd_stab_ltfu<- round(as.numeric(final.results$t0[[15]]),2)
mean_sd_stab_ltfu<- paste0(mean_stab_ltfu,"(",sd_stab_ltfu,")")
min_stab_ltfu<- round(as.numeric(final.results$t0[[16]]),2)
max_stab_ltfu<- round(as.numeric(final.results$t0[[17]]),2)
min_max_stab_ltfu<- paste0(min_stab_ltfu,"(",max_stab_ltfu,")")
#
#
#
# spec1 results (20 samples)---------------------------------------------------
spec1_ltfu_stab<- data.frame(
  Method=  c("HR","dHR","cHR","RR","cHR_km","RR_km"),
  Estimate= c(HR,dHR,cHR,RR,cHR_km,RR_km),
  LL= c(0.4761,0.4520,0.4601,0.6055,0.4938,0.6150),
  UL= c(0.6530,0.6218,0.6259,0.7365,1.0010,1.0006)
)
spec1_ltfu_stab[,sapply(spec1_ltfu_stab, is.numeric)]<- round(spec1_ltfu_stab[,sapply(spec1_ltfu_stab, is.numeric)],2)
write.csv(spec1_ltfu_stab,file="spec1_ltfu_stab.csv",row.names = T)
#
#
#
stabwts_spec1_ltfu<- cbind(mean_sd_stab_ltfu,min_max_stab_ltfu)
colnames(stabwts_spec1_ltfu)<- c("Mean(SD)","Min(Max)")
write.csv(stabwts_spec1_ltfu,file="stabwts_spec1_ltfu.csv",row.names = T)
#
#
#
# For results with unstab weights rerun the whole analysis 
# and replace weights= stabwt1 with weights=wt1
spec1_ltfu_unstab<- data.frame(
  Method=  c("HR","dHR","cHR","RR","cHR_km","RR_km"),
  Estimate= c(HR,dHR,cHR,RR,cHR_km,RR_km),
  LL= c(0.2737,0.3156,0.3234,0.4744,0.4938,0.6150),
  UL= c(0.7465,0.7172,0.7216,0.8132,1.0010,1.0006)
)
spec1_ltfu_unstab[,sapply(spec1_ltfu_unstab, is.numeric)]<- round(spec1_ltfu_unstab[,sapply(spec1_ltfu_unstab, is.numeric)],2)
write.csv(spec1_ltfu_unstab,file="spec1_ltfu_unstab.csv",row.names = T)
#
#
#
unstabwts_spec1_ltfu<- cbind(mean_sd_unstab_ltfu,min_max_unstab_ltfu)
colnames(unstabwts_spec1_ltfu)<- c("Mean(SD)","Min(Max)")
write.csv(unstabwts_spec1_ltfu,file="unstabwts_spec1_ltfu.csv",row.names = T)
#
#
#
# Forest Plots ltfu cHR---------------------------
# Analysis for Figure 3 in the manuscript
# For demonstration the results entered here are from the analysis reported in the manuscript
# and cannot be reproduced by running this script
data_forest<- data.frame(
  Specification= c("1","2","3","4","5","6"),
  mean= c(0.84, 0.82, 0.83, 0.83,0.81, 0.84),
  lower=c(0.53, 0.52, 0.52, 0.52,0.50, 0.53),
  upper=c(1.26, 1.23, 1.24, 1.25,1.21, 1.27),
  cHR= c("0.84", "0.82", "0.83", "0.83","0.81", "0.84"),
  CI= c("0.53,1.26","0.52,1.23","0.52,1.24","0.52,1.25","0.50,1.21","0.53,1.27")
)

forest_ltfu_cHR<- data_forest |>
  forestplot(labeltext= c(Specification,cHR,CI),
             #title= "Forest plot: Cumulative Hazard ratio with 95% Confidence interval of Overall Survival adjusted for Lost to follow-up estimated using IPCW"
             clip= c(0.2,2),
             vertices=TRUE,
             xlog=T,
             boxsize=0.1,
             xlab= "cHR (95%CI)",
             graph.pos=4,
             text_size=2,
             spacing=0,
             txt_gp=fpTxtGp(cex=0.95), # size
             axes=gpar(cex=0.7))|>
  fp_set_style(box="darkblue",
               line="black",
               summary=NULL)|>
  fp_add_header(Specification= c("","Specification"),
                cHR= c("","cHR"),
                CI= c("","95% CI"))|>
  fp_set_zebra_style("#EFEFEF")
forest_ltfu_cHR


#-------------------------------------------------------------------------------





# Analysis 2: XO adjusted analysis----------------------------------------------

# Analysis process for Table 4 and 5 in manuscript
#
#
#
# Commenting the code:
# In lines 761 to 960 entitled "Spec 1 to Spec 8"
# We provide the full code for the 8 analyses with explanatory comments 
# uncomment and execute only one specification at a time
# For demonstration we uncomment spec 1 and comment specs 2 to 8
# All subsequent steps (weight calculation and the different outcome models) should be run with any spec
# In lines 1188 to 1225 and for demonstration we provide the code for the forest plot with the results from the real data

df_xotd<- df_ipcw_pp[complete.cases(df_ipcw_pp$xotd),]

analyses.boot<- function(data, indices){
  
  ## Step 1: Generated bootstrap data, dd
  
  if (identical(indices, 1:length(indices))){
    dd <- data
  } else {
    dd <- create_bootstrap_data(data)
  }
  
  dd$wt2 <- NULL
  dd$stabwt2<- NULL
  
  knots0 = Hmisc::rcspline.eval(dd$time, knots.only = T, pc=T)
  knots1 = Hmisc::rcspline.eval(dd$AGE, knots.only = T, pc=T)
  knots2 = Hmisc::rcspline.eval(dd$tltd_wins, knots.only = T, pc=T)
  knots3 = Hmisc::rcspline.eval(dd$progtime, knots.only = T, pc=T)
  
  ## Step 2: Fit the weight models based on dd and calculate the weights
  ## Important: uncomment and execute only one specification at a time
  
  ### Spec 1----------------------------------------------------------------------
  # 
  ### Numerator control for XO
  nmodel_xotd_0<- glm( xotd ~ time, data=dd[complete.cases(dd$xotd) & dd$plannedtrt==0,],
                       family =binomial()) # Predict P for control arm

  # Predict
  dd$pnumerator2_0<- predict(nmodel_xotd_0, newdata=dd, type='response') # Set Numerator for Experimental =1 by setting P to 0
  dd$pnumerator2_1<- 0


  ### Denominator control for XO
  dmodel_xotd_0<- glm(
    xotd ~ time
    + SMKHISGR1 + ecogtdGR2 + AGE + tltd_wins + progtime
    + STRATRAN +  LIASE + DIAGINITGR3 + MICNSFL
    + I(progtime^2),
    data=dd[complete.cases(dd$xotd) & dd$plannedtrt==0 & dd$progtd==1,],
    family =binomial())

  # Predict
  dd$pdenominator2_0<- predict(dmodel_xotd_0, newdata=dd, type='response') # Set Denominator for Control non-progressors to 1 by setting P to 0
  dd$pdenominator2_0[dd$plannedtrt==0 & dd$progtd==0]=0 # Set Denominator for Experimental to 1 by setting P to 0
  dd$pdenominator2_1<- 0

  ### Spec 2----------------------------------------------------------------------
 #  
 #  ### Numerator control for XO
 #  nmodel_xotd_0<- glm( xotd ~ time, data=dd[complete.cases(dd$xotd) & dd$plannedtrt==0,],
 #                      family =binomial()) # Predict P for control arm
 # 
 # # Predict
 # dd$pnumerator2_0<- predict(nmodel_xotd_0, newdata=dd, type='response') # Set Numerator for Experimental =1 by setting P to 0
 # dd$pnumerator2_1<- 0
 # 
 # 
 #  ### Denominator control for XO
 #   dmodel_xotd_0<- glm(
 #     xotd ~ time
 #     + AGE + tltd_wins + progtime + I(progtime^2)
 #     + STRATRAN +  DIAGINITGR3,
 #     data=dd[complete.cases(dd$xotd) & dd$plannedtrt==0 & dd$progtd==1,],
 #     family =binomial())
 #   
 #  # Predict
 #  dd$pdenominator2_0<- predict(dmodel_xotd_0, newdata=dd, type='response') # Set Denominator for Control non-progressors to 1 by setting P to 0
 #  dd$pdenominator2_0[dd$plannedtrt==0 & dd$progtd==0]=0 # Set Denominator for Experimental to 1 by setting P to 0
 #  dd$pdenominator2_1<- 0
  
  ### Spec 3----------------------------------------------------------------------

  # ### Numerator control for XO
  # nmodel_xotd_0<- glm( xotd ~ rcs(time,knots0), data=dd[complete.cases(dd$xotd) & dd$plannedtrt==0,],
  #                      family =binomial()) # Predict P for control arm
  # 
  # # Predict
  # dd$pnumerator2_0<- predict(nmodel_xotd_0, newdata=dd, type='response') # Set Numerator for Experimental =1 by setting P to 0
  # dd$pnumerator2_1<- 0
  # 
  # 
  # ### Denominator control for XO
  # dmodel_xotd_0<- glm(
  #   xotd ~ rcs(time,knots0)
  #   + SMKHISGR1 + ecogtdGR2
  #   + rcs(AGE,knots1)
  #   + rcs(tltd_wins,knots2)
  #   + rcs(progtime,knots3)
  #   + STRATRAN +  LIASE + DIAGINITGR3 + MICNSFL,
  #   data=dd[complete.cases(dd$xotd) & dd$plannedtrt==0 & dd$progtd==1,],
  #   family =binomial())
  # 
  # # Predict
  # dd$pdenominator2_0<- predict(dmodel_xotd_0, newdata=dd, type='response') # Set Denominator for Control non-progressors to 1 by setting P to 0
  # dd$pdenominator2_0[dd$plannedtrt==0 & dd$progtd==0]=0 # Set Denominator for Experimental to 1 by setting P to 0
  # dd$pdenominator2_1<- 0
  # 
  ### Spec 4----------------------------------------------------------------------
  # # 
  # ### Numerator control for XO
  # nmodel_xotd_0<- glm( xotd ~ time, data=dd[complete.cases(dd$xotd) & dd$plannedtrt==0,],
  #                      family =binomial()) # Predict P for control arm
  # 
  # # Predict
  # dd$pnumerator2_0<- predict(nmodel_xotd_0, newdata=dd, type='response') # Set Numerator for Experimental =1 by setting P to 0
  # dd$pnumerator2_1<- 0
  # 
  # 
  # ### Denominator control for XO
  # dmodel_xotd_0<- glm(
  #   xotd ~ time
  #   + SMKHISGR1 + ecogtdGR1 + AGE + tltd_wins + progtime
  #   + STRATRAN +  LIASE + DIAGINITGR3 + MICNSFL
  #   + I(progtime^2),
  #   data=dd[complete.cases(dd$xotd) & dd$plannedtrt==0 & dd$progtd==1,],
  #   family =binomial())
  # 
  # # Predict
  # dd$pdenominator2_0<- predict(dmodel_xotd_0, newdata=dd, type='response') # Set Denominator for Control non-progressors to 1 by setting P to 0
  # dd$pdenominator2_0[dd$plannedtrt==0 & dd$progtd==0]=0 # Set Denominator for Experimental to 1 by setting P to 0
  # dd$pdenominator2_1<- 0
  
  ### Spec 5----------------------------------------------------------------------
  
  # ### Numerator control for XO
  # nmodel_xotd_0<- glm( xotd ~ rcs(time,knots0), data=dd[complete.cases(dd$xotd) & dd$plannedtrt==0,],
  #                      family =binomial()) # Predict P for control arm
  # 
  # # Predict
  # dd$pnumerator2_0<- predict(nmodel_xotd_0, newdata=dd, type='response') # Set Numerator for Experimental =1 by setting P to 0
  # dd$pnumerator2_1<- 0
  # 
  # 
  # ### Denominator control for XO
  # dmodel_xotd_0<- glm(
  #   xotd ~ rcs(time,knots0)
  #   + rcs(AGE,knots1)
  #   + rcs(tltd_wins,knots2)
  #   + rcs(progtime,knots3)
  #   + STRATRAN +  DIAGINITGR3,
  #   data=dd[complete.cases(dd$xotd) & dd$plannedtrt==0 & dd$progtd==1,],
  #   family =binomial())
  # 
  # # Predict
  # dd$pdenominator2_0<- predict(dmodel_xotd_0, newdata=dd, type='response') # Set Denominator for Control non-progressors to 1 by setting P to 0
  # dd$pdenominator2_0[dd$plannedtrt==0 & dd$progtd==0]=0 # Set Denominator for Experimental to 1 by setting P to 0
  # dd$pdenominator2_1<- 0
  
  ### Spec 6----------------------------------------------------------------------

  # ### Numerator control for XO
  # nmodel_xotd_0<- glm( xotd ~ time, data=dd[complete.cases(dd$xotd) & dd$plannedtrt==0,],
  #                      family =binomial()) # Predict P for control arm
  # 
  # # Predict
  # dd$pnumerator2_0<- predict(nmodel_xotd_0, newdata=dd, type='response') # Set Numerator for Experimental =1 by setting P to 0
  # dd$pnumerator2_1<- 0
  # 
  # 
  # ### Denominator control for XO
  # dmodel_xotd_0<- glm(
  #   xotd ~ time
  #   + AGE + tltd_wins + progtime + I(progtime^2)
  #   + STRATRAN,
  #   data=dd[complete.cases(dd$xotd) & dd$plannedtrt==0 & dd$progtd==1,],
  #   family =binomial())
  # 
  # # Predict
  # dd$pdenominator2_0<- predict(dmodel_xotd_0, newdata=dd, type='response') # Set Denominator for Control non-progressors to 1 by setting P to 0
  # dd$pdenominator2_0[dd$plannedtrt==0 & dd$progtd==0]=0 # Set Denominator for Experimental to 1 by setting P to 0
  # dd$pdenominator2_1<- 0
  
  ### Spec 7----------------------------------------------------------------------
  
  # ### Numerator control for XO
  # nmodel_xotd_0<- glm( xotd ~ time, data=dd[complete.cases(dd$xotd) & dd$plannedtrt==0,],
  #                      family =binomial()) # Predict P for control arm
  # 
  # # Predict
  # dd$pnumerator2_0<- predict(nmodel_xotd_0, newdata=dd, type='response') # Set Numerator for Experimental =1 by setting P to 0
  # dd$pnumerator2_1<- 0
  # 
  # 
  # ### Denominator control for XO
  # dmodel_xotd_0<- glm(
  #   xotd ~ time
  #   + SMKHISGR1 + ecogtdGR2 + AGE + tltd_wins
  #   + STRATRAN +  LIASE + DIAGINITGR3 + MICNSFL,
  #   data=dd[complete.cases(dd$xotd) & dd$plannedtrt==0 & dd$progtd==1,],
  #   family =binomial())
  # 
  # # Predict
  # dd$pdenominator2_0<- predict(dmodel_xotd_0, newdata=dd, type='response') # Set Denominator for Control non-progressors to 1 by setting P to 0
  # dd$pdenominator2_0[dd$plannedtrt==0 & dd$progtd==0]=0 # Set Denominator for Experimental to 1 by setting P to 0
  # dd$pdenominator2_1<- 0
  
  ### Spec 8----------------------------------------------------------------------
  
  # ### Numerator control for XO
  # nmodel_xotd_0<- glm( xotd ~ time, data=dd[complete.cases(dd$xotd) & dd$plannedtrt==0,],
  #                      family =binomial()) # Predict P for control arm
  # 
  # # Predict
  # dd$pnumerator2_0<- predict(nmodel_xotd_0, newdata=dd, type='response') # Set Numerator for Experimental =1 by setting P to 0
  # dd$pnumerator2_1<- 0
  # 
  # 
  # ### Denominator control for XO
  # dmodel_xotd_0<- glm(
  #   xotd ~ time
  #   + AGE + tltd_wins
  #   + STRATRAN +  DIAGINITGR3,
  #   data=dd[complete.cases(dd$xotd) & dd$plannedtrt==0 & dd$progtd==1,],
  #   family =binomial())
  # 
  # # Predict
  # dd$pdenominator2_0<- predict(dmodel_xotd_0, newdata=dd, type='response') # Set Denominator for Control non-progressors to 1 by setting P to 0
  # dd$pdenominator2_0[dd$plannedtrt==0 & dd$progtd==0]=0 # Set Denominator for Experimental to 1 by setting P to 0
  # dd$pdenominator2_1<- 0
  
  ### Weight Calculation----------------------------------------------------------
  
  dd<- dd[order(dd$id,dd$time),] # Contribution from xotd(t)=0 is 1-P # Contribution from xotd(t)=1 is P
    
  dd<- dd%>%
    dplyr::mutate(num2= ifelse(plannedtrt==1,
                               (xotd*pnumerator2_1 + (1-xotd)*(1-pnumerator2_1)),
                               (xotd*pnumerator2_0 + (1-xotd)* (1-pnumerator2_0))),        
                  den2=ifelse(plannedtrt==1,
                              (xotd*pdenominator2_1 + (1-xotd)* (1-pdenominator2_1)),
                              (xotd*pdenominator2_0 + (1-xotd)* (1-pdenominator2_0))),
                  num2= ifelse(time==0,1,num2),
                  den2= ifelse(time==0,1,den2)
    ) %>%
    dplyr:: group_by(id) %>%
    dplyr:: mutate(
      n2= cumprod(num2),
      d2= cumprod(den2)
    )%>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      stabwt2= n2/d2,
      wt2=1/d2
    )
  
  dd_wts_xo= dd[complete.cases(dd$xotd) & dd$xotd==0,]
  
  # Distribution of unstab wts
  mean_unstab_xo<- mean(dd_wts_xo$wt2)
  sd_unstab_xo<- sd(dd_wts_xo$wt2)
  min_unstab_xo<- min(dd_wts_xo$wt2)
  max_unstab_xo<- max(dd_wts_xo$wt2)
  
  
  # Distribution of stab wts
  mean_stab_xo<- mean(dd_wts_xo$stabwt2)
  sd_stab_xo<- sd(dd_wts_xo$stabwt2)
  min_stab_xo<- min(dd_wts_xo$stabwt2)
  max_stab_xo<- max(dd_wts_xo$stabwt2)
  
  ## Step 3: Remove observations at censoring and fit the outcomes model--------
  
  dd_outcome<- dd[dd$xotd==0,]
  
  # Cox proportional hazard model
  
  xo.cox <- coxph(Surv(time,time+1,deathtd) ~ plannedtrt + cluster(id), robust = T,
                    data = dd_outcome, weights= stabwt2, ties = ("breslow")) # weights=wt2
  
  ### HR
  HR<- exp(coef(xo.cox))[1]
  
  
  # Pooled logistic regression (dHR, cHR, RR)
  xo.logit<- glm(deathtd ~ time + time2 + factor(plannedtrt),
                   data= dd_outcome, weights= stabwt2, family =quasibinomial()) #weights=wt2
  
  ### dHR
  dHR<- exp(coef(xo.logit))[4]
  ### cHR and RR
  dd_outcome$prob<- 1-predict(xo.logit, newdata=dd_outcome, type="response") # prob: survival probability
  dd_1<- dd_outcome%>%
    dplyr::arrange(id,time)%>%
    dplyr::group_by(id)%>%
    dplyr::mutate(
      surv=cumprod(prob))
  dd_1<- dd_1[,c('surv','plannedtrt','time')]
  results<- dd_1%>%
    dplyr::group_by(time, plannedtrt)%>%
    dplyr::summarize(mean_survival=mean(surv))
  results<- results%>%
    dplyr::ungroup()%>%
    dplyr:: mutate(
      time=time+1
    )
  results<- dplyr::bind_rows(c(time=0,plannedtrt=0,mean_survival=1),c(time=0,plannedtrt=1,mean_survival=1),results )
  results$plannedtrtf<- factor(results$plannedtrt,labels=c("Crizotinib","Brigatinib"))
  wide_results<- reshape2::dcast(results, time ~ plannedtrtf,value.var = 'mean_survival')
  wide_results<- wide_results%>%
    mutate(
      risk0=(1-Crizotinib),
      risk1=(1-Brigatinib),
      RD= (1-Brigatinib)-(1-Crizotinib),
      cHR= log(Brigatinib)/log(Crizotinib),
      RR=(1-Brigatinib)/(1-Crizotinib))
  
  
  # KM Estimator (Note: run K= 36)
  kmfit<- survfit(Surv(time, time+1,deathtd) ~ plannedtrt, weights = stabwt2, data = dd_outcome) # weights=wt2
  
  ### cHR
  S0<- summary(kmfit, times=k, extend=T)$surv[1]
  S1<- summary(kmfit, times=k, extend=T)$surv[2]
  cHR_km<- log(S1)/log(S0)
  
  ### RR
  R0_km<- 1-summary(kmfit, times = k)$surv[1]
  R1_km<- 1-summary(kmfit, times = k)$surv[2]
  RD_km<- R1_km-R0_km
  RR_km<- R1_km/R0_km

  
  # Return results--------------------------------------------------------------
  
  #assign("data_xo_final",dd_outcome,envir = .GlobalEnv)
  #assign("results_xo_final",results,envir = .GlobalEnv)
  
  return(c(HR,
           dHR,
           wide_results$risk0[which(wide_results$time==k-1)],
           wide_results$risk1[which(wide_results$time==k-1)],
           wide_results$RD[which(wide_results$time==k-1)],
           wide_results$cHR[which(wide_results$time==k-1)],
           wide_results$RR[which(wide_results$time==k-1)],
           cHR_km,
           RR_km,
           mean_unstab_xo,
           sd_unstab_xo,
           min_unstab_xo,
           max_unstab_xo,
           mean_stab_xo,
           sd_stab_xo,
           min_stab_xo,
           max_stab_xo))
} 

# Set replicate failures to NA
boot.robust<- function(data,indices){
  tryCatch({analyses.boot(data,indices)}, error= function(e){NA})
}

set.seed(7268)

final.results<- boot(data= df_xotd, statistic= boot.robust, R=20) 
final.results$t0

# CoxPH
HR<- round(final.results$t0[[1]],2)
HR
conf_interval_HR <- boot.ci(final.results, conf = 0.95, type = "perc", index = 1)
conf_interval_HR


# Pooled logistic regression
dHR<- round(final.results$t0[[2]],2)
dHR
conf_interval_dHR <- boot.ci(final.results, conf = 0.95, type = "perc", index = 2)
conf_interval_dHR
#
cHR<- round(final.results$t0[[6]],2)
cHR
conf_interval_cHR <- boot.ci(final.results, conf = 0.95, type = "perc", index = 6)
conf_interval_cHR
#
RR<- round(final.results$t0[[7]],2)
RR
conf_interval_RR <- boot.ci(final.results, conf = 0.95, type = "perc", index = 7)
conf_interval_RR


# KM Estimator
cHR_km<-round(final.results$t0[[8]],2)
cHR_km
conf_interval_cHR_km<- boot.ci(final.results, conf = 0.95, type = "perc", index = 8)
conf_interval_cHR_km
#
RR_km<-round(final.results$t0[[9]],2)
RR_km
conf_interval_RR_km<- boot.ci(final.results, conf = 0.95, type = "perc", index = 9)
conf_interval_RR_km
#
#
#
# Distribution of XO unstab wts
mean_unstab_xo<- round(as.numeric(final.results$t0[[10]]),2)
sd_unstab_xo<- round(as.numeric(final.results$t0[[11]]),2)
mean_sd_unstab_xo<- paste0(mean_unstab_xo,"(",sd_unstab_xo,")")
min_unstab_xo<- round(final.results$t0[[12]],2)
max_unstab_xo<- round(final.results$t0[[13]],2)
min_max_unstab_xo<- paste0(min_unstab_xo,"(",max_unstab_xo,")")

# Distribution of xo stab wts
mean_stab_xo<- round(as.numeric(final.results$t0[[14]]),2)
sd_stab_xo<- round(as.numeric(final.results$t0[[15]]),2)
mean_sd_stab_xo<- paste0(mean_stab_xo,"(",sd_stab_xo,")")
min_stab_xo<- round(as.numeric(final.results$t0[[16]]),2)
max_stab_xo<- round(as.numeric(final.results$t0[[17]]),2)
min_max_stab_xo<- paste0(min_stab_xo,"(",max_stab_xo,")")
#
#
#
# spec1 results (20 samples)--------------------------
spec1_xo_stab<- data.frame(
  Method=  c("HR","dHR","cHR","RR","cHR_km","RR_km"),
  Estimate= c(HR,dHR,cHR,RR,cHR_km,RR_km),
  LL= c(0.5091,0.4609,0.4688,0.607,0.3377,0.5000),
  UL= c(1.9743,2.0691,2.0577,1.828,1.2784,1.2027 )
)
spec1_xo_stab[,sapply(spec1_xo_stab, is.numeric)]<- round(spec1_xo_stab[,sapply(spec1_xo_stab, is.numeric)],2)
write.csv(spec1_xo_stab,file="spec1_xo_stab.csv",row.names = T)
#
#
#
stabwts_spec1_xo<- cbind(mean_sd_stab_xo,min_max_stab_xo)
colnames(stabwts_spec1_xo)<- c("Mean(SD)","Min(Max)")
write.csv(stabwts_spec1_xo,file="stabwts_spec1_xo.csv",row.names = T)
#
#
#
# For results with unstab weights rerun the whole analysis 
# and replace weights= stabwt2 with weights=wt2
spec1_xo_unstab<- data.frame(
  Method=  c("HR","dHR","cHR","RR","cHR_km","RR_km"),
  Estimate= c(HR,dHR,cHR,RR,cHR_km,RR_km),
  LL= c(0.5331,0.4614,0.4712,0.6181,0.3377,0.5000),
  UL= c(1.8496,1.9681,1.9566,1.7406,1.2784,1.2027 )
)
spec1_xo_unstab[,sapply(spec1_xo_unstab, is.numeric)]<- round(spec1_xo_unstab[,sapply(spec1_xo_unstab, is.numeric)],2)
write.csv(spec1_xo_unstab,file="spec1_xo_unstab.csv",row.names = T)
#
#
#
unstabwts_spec1_xo<- cbind(mean_sd_unstab_xo,min_max_unstab_xo)
colnames(unstabwts_spec1_xo)<- c("Mean(SD)","Min(Max)")
write.csv(unstabwts_spec1_xo,file="unstabwts_spec1_xo.csv",row.names = T)



# Forest Plots xo cHR---------------------------
# Analysis for Figure 4 in the manuscript
# For demonstration the results entered here are from the analysis reported in the manuscript
# and cannot be reproduced by running this script
data_forest<- data.frame(
  Specification= c("1","2","3","4","5","6","7","8"),
  mean= c(0.69,0.69, 0.67,0.69,0.65,0.73,0.40,0.38),
  lower=c(0.39,0.41,0.38,0.39,0.32,0.45,0.11,0.12),
  upper=c(1.21, 1.21, 1.26,1.21,1.22,1.22, 1.20, 0.98),
  cHR= c("0.69","0.69", "0.67" ,"0.69","0.65","0.73","0.40","0.38"),
  CI= c("0.39,1.21","0.41,1.21","0.38,1.26","0.39,1.21","0.32,1.22","0.45,1.22",
        "0.11,1.20","0.12,0.98")
)

forest_xo_cHR<- data_forest |>
  forestplot(labeltext= c(Specification,cHR,CI),
             #title= "Forest plot: Cumulative Hazard ratio with 95% Confidence interval of Overall Survival adjusted for Lost to follow-up estimated using IPCW"
             clip= c(0.1,2),
             vertices=TRUE,
             xlog=T,
             boxsize=0.1,
             xlab= "cHR (95%CI)",
             graph.pos=4,
             text_size=2,
             spacing=0,
             txt_gp=fpTxtGp(cex=0.95), # size
             axes=gpar(cex=0.7))|>
  fp_set_style(box="darkblue",
               line="black",
               summary=NULL)|>
  fp_add_header(Specification= c("","Specification"),
                cHR= c("","cHR"),
                CI= c("","95% CI"))|>
  fp_set_zebra_style("#EFEFEF")
forest_xo_cHR



# XO adjusted analysis "Spec 7" with different truncation --------------------
# Analysis for Table 5 in the manuscript 
# Note: uncomment and execute only one truncation percentile at a time 
# here we uncomment truncation percentile: 95

df_xotd<- df_ipcw_pp[complete.cases(df_ipcw_pp$xotd),]

analyses.boot<- function(data, indices){
  
  ## Step 1: Generated bootstrap data, dd
  if (identical(indices, 1:length(indices))){
    dd <- data
  } else {
    dd <- create_bootstrap_data(data)
  }
  
  dd$wt2 <- NULL
  dd$stabwt2<- NULL
  
  
  ## Step 2: Fit the weight models based on dd and calculate the weights
  
  ### Numerator control for XO
  nmodel_xotd_0<- glm( xotd ~ time, data=dd[complete.cases(dd$xotd) & dd$plannedtrt==0,],
                       family =binomial()) # Predict P for control arm

  # Predict
  dd$pnumerator2_0<- predict(nmodel_xotd_0, newdata=dd, type='response') # Set Numerator for Experimental =1 by setting P to 0
  dd$pnumerator2_1<- 0


  ### Denominator control for XO
  dmodel_xotd_0<- glm(
    xotd ~ time
    + SMKHISGR1 + ecogtdGR2 + AGE + tltd_wins
    + STRATRAN +  LIASE + DIAGINITGR3 + MICNSFL,
    data=dd[complete.cases(dd$xotd) & dd$plannedtrt==0 & dd$progtd==1,],
    family =binomial())

  # Predict
  dd$pdenominator2_0<- predict(dmodel_xotd_0, newdata=dd, type='response') # Set Denominator for Control non-progressors to 1 by setting P to 0
  dd$pdenominator2_0[dd$plannedtrt==0 & dd$progtd==0]=0 # Set Denominator for Experimental to 1 by setting P to 0
  dd$pdenominator2_1<- 0
  

  ### Weight Calculation
  
  dd<- dd[order(dd$id,dd$time),] # Contribution from xotd(t)=0 is 1-P # Contribution from xotd(t)=1 is P
  
  dd<- dd%>%
    dplyr::mutate(num2= ifelse(plannedtrt==1,
                               (xotd*pnumerator2_1 + (1-xotd)*(1-pnumerator2_1)),
                               (xotd*pnumerator2_0 + (1-xotd)* (1-pnumerator2_0))),        
                  den2=ifelse(plannedtrt==1,
                              (xotd*pdenominator2_1 + (1-xotd)* (1-pdenominator2_1)),
                              (xotd*pdenominator2_0 + (1-xotd)* (1-pdenominator2_0))),
                  num2= ifelse(time==0,1,num2),
                  den2= ifelse(time==0,1,den2)
    ) %>%
    dplyr:: group_by(id) %>%
    dplyr:: mutate(
      n2= cumprod(num2),
      d2= cumprod(den2)
    )%>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      stabwt2= n2/d2,
      wt2=1/d2
    )
  
  # # truncation percentile: 100
  # percentile_100 <- quantile(dd$stabwt2, 1)
  # dd$tstabwt2<- dd$stabwt2
  # dd$tstabwt2[dd$stabwt2>percentile_100]<- percentile_100
  # 
  # # truncation percentile: 99
  # percentile_99 <- quantile(dd$stabwt2, 0.99)
  # dd$tstabwt2<- dd$stabwt2
  # dd$tstabwt2[dd$stabwt2>percentile_99]<- percentile_99
  # 
  # truncation percentile: 95
  percentile_95 <- quantile(dd$stabwt2, 0.95)
  dd$tstabwt2<- dd$stabwt2
  dd$tstabwt2[dd$stabwt2>percentile_95]<- percentile_95
  # 
  # # truncation percentile: 90
  # percentile_90 <- quantile(dd$stabwt2, 0.90)
  # dd$tstabwt2<- dd$stabwt2
  # dd$tstabwt2[dd$stabwt2>percentile_90]<- percentile_90
  # # 
  # # truncation percentile: 75
  # percentile_75 <- quantile(dd$stabwt2, 0.75)
  # dd$tstabwt2<- dd$stabwt2
  # dd$tstabwt2[dd$stabwt2>percentile_75]<- percentile_75
  # 
  # # truncation percentile: 50
  # percentile_50 <- quantile(dd$stabwt2, 0.5)
  # dd$tstabwt2<- dd$stabwt2
  # dd$tstabwt2[dd$stabwt2>percentile_50]<- percentile_50

  # Distribution of truncated weights
  data_all= dd
  data= data_all[data_all$xotd==0,]
  mean_tstab<- round(mean(data$tstabwt2),2)
  sd_tstab<- round(sd(data$tstabwt2),2)
  min_tstab<- round(min(data$tstabwt2),2)
  max_tstab<- round(max(data$tstabwt2),2)
  
  ## Step 3: Remove observations at censoring and fit the outcomes model
  dd_outcome<- dd[dd$xotd==0,]
  
  # Cox proportional hazard model
  
  xo.cox <- coxph(Surv(time,time+1,deathtd) ~ plannedtrt + cluster(id), robust = T,
                  data = dd_outcome, weights=tstabwt2, ties = ("breslow"))
  
  ### HR
  HR<- exp(coef(xo.cox))[1]
  
  
  # Pooled logistic regression (dHR, cHR, RR)
  xo.logit<- glm(deathtd ~ time + time2 + factor(plannedtrt),
                 data= dd_outcome, weights= tstabwt2, family =quasibinomial())
  
  ### dHR
  dHR<- exp(coef(xo.logit))[4]
  ### cHR and RR
  dd_outcome$prob<- 1-predict(xo.logit, newdata=dd_outcome, type="response") # prob: survival probability
  dd_1<- dd_outcome%>%
    dplyr::arrange(id,time)%>%
    dplyr::group_by(id)%>%
    dplyr::mutate(
      surv=cumprod(prob))
  dd_1<- dd_1[,c('surv','plannedtrt','time')]
  results<- dd_1%>%
    dplyr::group_by(time, plannedtrt)%>%
    dplyr::summarize(mean_survival=mean(surv))
  results<- results%>%
    dplyr::ungroup()%>%
    dplyr:: mutate(
      time=time+1
    )
  results<- dplyr::bind_rows(c(time=0,plannedtrt=0,mean_survival=1),c(time=0,plannedtrt=1,mean_survival=1),results )
  results$plannedtrtf<- factor(results$plannedtrt,labels=c("Crizotinib","Brigatinib"))
  wide_results<- reshape2::dcast(results, time ~ plannedtrtf,value.var = 'mean_survival')
  wide_results<- wide_results%>%
    mutate(
      risk0=(1-Crizotinib),
      risk1=(1-Brigatinib),
      RD= (1-Brigatinib)-(1-Crizotinib),
      cHR= log(Brigatinib)/log(Crizotinib),
      RR=(1-Brigatinib)/(1-Crizotinib))
  
  
  # KM Estimator (Note: run K= 36)
  kmfit<- survfit(Surv(time, time+1,deathtd) ~ plannedtrt, weights = tstabwt2, data = dd_outcome)
  
  ### cHR
  S0<- summary(kmfit, times=k, extend=T)$surv[1]
  S1<- summary(kmfit, times=k, extend=T)$surv[2]
  cHR_km<- log(S1)/log(S0)
  
  ### RR
  R0_km<- 1-summary(kmfit, times = k)$surv[1]
  R1_km<- 1-summary(kmfit, times = k)$surv[2]
  RD_km<- R1_km-R0_km
  RR_km<- R1_km/R0_km
  
  
  # Return results
  #print(head(dd_outcome))
  assign("data_xo_final",dd_outcome,envir = .GlobalEnv)
  assign("results_xo_final",results,envir = .GlobalEnv)
  
  return(c(HR,
           dHR,
           wide_results$risk0[which(wide_results$time==k-1)],
           wide_results$risk1[which(wide_results$time==k-1)],
           wide_results$RD[which(wide_results$time==k-1)],
           wide_results$cHR[which(wide_results$time==k-1)],
           wide_results$RR[which(wide_results$time==k-1)],
           cHR_km,
           RR_km,
           mean_tstab,
           sd_tstab,
           min_tstab,
           max_tstab))
} 


boot.robust<- function(data,indices){
  tryCatch({analyses.boot(data,indices)}, error= function(e){NA})
}

set.seed(7268)
final.results<- boot(data= df_xotd, statistic= boot.robust, R=20) 
final.results$t0

# CoxPH
HR<- round(final.results$t0[[1]],2)
HR
conf_interval_HR <- boot.ci(final.results, conf = 0.95, type = "perc", index = 1)
conf_interval_HR


# Pooled logistic regression
dHR<- round(final.results$t0[[2]],2)
dHR
conf_interval_dHR <- boot.ci(final.results, conf = 0.95, type = "perc", index = 2)
conf_interval_dHR
#
cHR<- round(final.results$t0[[6]],2)
cHR
conf_interval_cHR <- boot.ci(final.results, conf = 0.95, type = "perc", index = 6)
conf_interval_cHR
#
RR<- round(final.results$t0[[7]],2)
RR
conf_interval_RR <- boot.ci(final.results, conf = 0.95, type = "perc", index = 7)
conf_interval_RR


# KM Estimator
cHR_km<-round(final.results$t0[[8]],2)
cHR_km
conf_interval_cHR_km<- boot.ci(final.results, conf = 0.95, type = "perc", index = 8)
conf_interval_cHR_km
#
RR_km<-round(final.results$t0[[9]],2)
RR_km
conf_interval_RR_km<- boot.ci(final.results, conf = 0.95, type = "perc", index = 9)
conf_interval_RR_km
#
# Distribution of xo tstab wts
mean_tstab<- round(as.numeric(final.results$t0[[10]]),2)
sd_tstab<- round(as.numeric(final.results$t0[[11]]),2)
mean_sd_tstab<- paste0(mean_tstab,"(",sd_tstab,")")
min_tstab<- round(as.numeric(final.results$t0[[12]]),2)
max_tstab<- round(as.numeric(final.results$t0[[13]]),2)
min_max_tstab<- paste0(min_tstab,"(",max_tstab,")")
#
#
# results for truncation percentile 95 (20 samples) ---------------------------
spec7_tstab_95<- data.frame(
  Method=  c("HR","dHR","cHR","RR","cHR_km","RR_km"),
  Estimate= c(HR,dHR,cHR,RR,cHR_km,RR_km),
  LL= c(1.254,1.260,1.258,1.205,1.082,1.064),
  UL= c(1.881,1.898,1.886,1.660,1.636,1.466)
)
spec7_tstab_95[,sapply(spec7_tstab_95, is.numeric)]<- round(spec7_tstab_95[,sapply(spec7_tstab_95, is.numeric)],2)
write.csv(spec7_tstab_95,file="spec7_tstab_95.csv",row.names = T)
#
tstabwt_spec7_95<- cbind(mean_sd_tstab,min_max_tstab)
colnames(tstabwt_spec7_95)<- c("Mean(SD)","Min(Max)")
write.csv(tstabwt_spec7_95,file="tstabwt_spec7_95.csv",row.names = T)
#
#
#-------------------------------------------------------------------------------

# Analysis 3: LTFU and XO adjusted analysis-------------------------------------

# Analysis process for Table 3 IPCW approach in the manuscript
### Important: uncomment and execute only one specification at a time
### For demonstration we uncomment spec 4 and comment spec 7

df_censoring_xotd <- df_ipcw_pp

analyses.boot<- function(data, indices){
  
  ## Step 1: Generated bootstrap data, dd
  
  if (identical(indices, 1:length(indices))){
    dd <- data
  } else {
    dd <- create_bootstrap_data(data)
  }
  
  dd$wt1 <- NULL
  dd$stabwt1<- NULL
  dd$wt2 <- NULL
  dd$stabwt2<- NULL
  
  ## Step 2: Weighting model for LTFU and calculate the weights
  
  ### Numerator control for LTFU
  nmodel_lofu_0<- glm(censoringtd_new ~ time + time2,
                      data=dd[dd$plannedtrt==0,], family =binomial())
  
  ### Numerator exp for LTFU
  nmodel_lofu_1<- glm(censoringtd_new ~ time + time2, 
                      data=dd[dd$plannedtrt==1,], family =binomial())
  
  # Predict
  dd$pnumerator1_0<- predict(nmodel_lofu_0, newdata=dd, type='response')
  dd$pnumerator1_1<- predict(nmodel_lofu_1, newdata=dd, type='response')
  
  ### Denominator control for LTFU
  dmodel_lofu_0<- glm(
    censoringtd_new ~ time + SMKHISGR1 + STRATRAN + ecogtdGR1
    + AGE + MICNSFL + tltd_wins + progtd + trttd + LIASE + DIAGINITGR3
    + time2
    + AGE2,
    data=dd[dd$plannedtrt==0,],
    family =binomial())
  
  ### Denominator exp for LTFU
  dmodel_lofu_1<- glm(
    censoringtd_new ~ time + AGE + SEX + LIASE + STRATRAN + ECOGGR1 + ecogtdGR1 
    + tltd_wins + icprogtd + progtd + RACEGR1 + SMKHISGR1 + MICNSFL 
    + DIAGINITGR3 + PRRADYN + time2 + AGE2,
    data=dd[dd$plannedtrt==1,],
    family =binomial())
  
  # Predict
  dd$pdenominator1_0<- predict(dmodel_lofu_0, newdata=dd, type='response')
  dd$pdenominator1_1<- predict(dmodel_lofu_1, newdata=dd, type='response')
  
  
  ### Weight Calculation for LTFU
  
  # Order the data by "id" and "time"
  dd<- dd[order(dd$id,dd$time),]
  
  # Calculate num1 and den1 
  dd<- dd%>%dplyr::mutate(
    num1= ifelse(plannedtrt==1,
                 (censoringtd_new*pnumerator1_1 + (1-censoringtd_new)*(1-pnumerator1_1)),
                 (censoringtd_new*pnumerator1_0 + (1-censoringtd_new)* (1-pnumerator1_0))),
    den1=ifelse(plannedtrt==1,
                (censoringtd_new*pdenominator1_1 + (1-censoringtd_new)* (1-pdenominator1_1)),
                (censoringtd_new*pdenominator1_0 + (1-censoringtd_new)* (1-pdenominator1_0))),
    num1= ifelse(time==0,1,num1),
    den1= ifelse(time==0,1,den1)
  ) %>%
    dplyr::group_by(id) %>%
    dplyr::mutate(
      n1= cumprod(num1),
      d1= cumprod(den1)
    )%>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      stabwt1= n1/d1, # stabilized weights
      wt1=1/d1        # unstabilized weights 
    )
  
  dd_wts= dd[dd$censoringtd_new==0,]
  
  # Distribution of ltfu unstab wts
  mean_unstab<- mean(dd_wts$wt1)
  sd_unstab<- sd(dd_wts$wt1)
  min_unstab<- min(dd_wts$wt1)
  max_unstab<- max(dd_wts$wt1)
  
  
  # Distribution of ltfu Stabwts
  mean_stab<- mean(dd_wts$stabwt1)
  sd_stab<- sd(dd_wts$stabwt1)
  min_stab<- min(dd_wts$stabwt1)
  max_stab<- max(dd_wts$stabwt1)
  
  ## Step 3: Weighting model for LTFU
  ## Note: uncomment and execute only one spec at a time (either 4 or 7)
  
  ### spec 4 --------------
  
  ### Numerator control for XO
  nmodel_xotd_0<- glm( xotd ~ time, data=dd[complete.cases(dd$xotd) & dd$plannedtrt==0,], family =binomial()) # Predict P for control arm

  # Predict
  dd$pnumerator2_0<- predict(nmodel_xotd_0, newdata=dd, type='response') #Set Numerator for Experimental =1 by setting P to 0
  dd$pnumerator2_1<- 0


  ### Denominator control for XO
  dmodel_xotd_0<- glm(
    xotd ~ time
    + SMKHISGR1 + ecogtdGR1 + AGE + tltd_wins + progtime
    + STRATRAN +  LIASE + DIAGINITGR3 + MICNSFL
    + I(progtime^2),
    data=dd[complete.cases(dd$xotd) & dd$plannedtrt==0 & dd$progtd==1,],
    family =binomial())

  # Predict
  dd$pdenominator2_0<- predict(dmodel_xotd_0, newdata=dd, type='response') # Set Denominator for Control non-progressors to 1 by setting P to 0
  dd$pdenominator2_0[dd$plannedtrt==0 & dd$progtd==0]=0 # Set Denominator for Experimental to 1 by setting P to 0
  dd$pdenominator2_1<- 0
  
  ### spec 7 --------------
  
  # ### Numerator control for XO
  # nmodel_xotd_0<- glm( xotd ~ time, data=dd[complete.cases(dd$xotd) & dd$plannedtrt==0,],
  #                      family =binomial()) # Predict P for control arm
  # 
  # # Predict
  # dd$pnumerator2_0<- predict(nmodel_xotd_0, newdata=dd, type='response') # Set Numerator for Experimental =1 by setting P to 0
  # dd$pnumerator2_1<- 0
  # 
  # 
  # ### Denominator control for XO
  # dmodel_xotd_0<- glm(
  #   xotd ~ time
  #   + SMKHISGR1 + ecogtdGR2 + AGE + tltd_wins
  #   + STRATRAN +  LIASE + DIAGINITGR3 + MICNSFL,
  #   data=dd[complete.cases(dd$xotd) & dd$plannedtrt==0 & dd$progtd==1,],
  #   family =binomial())
  # 
  # # Predict
  # dd$pdenominator2_0<- predict(dmodel_xotd_0, newdata=dd, type='response') # Set Denominator for Control non-progressors to 1 by setting P to 0
  # dd$pdenominator2_0[dd$plannedtrt==0 & dd$progtd==0]=0 # Set Denominator for Experimental to 1 by setting P to 0
  # dd$pdenominator2_1<- 0
  # 
  
  # Weight Calculation for XO
  
  dd<- dd[order(dd$id,dd$time),] 
  
  dd<- dd%>%
    dplyr::mutate(num2= ifelse(plannedtrt==1,
                               (xotd*pnumerator2_1 + (1-xotd)*(1-pnumerator2_1)),
                               (xotd*pnumerator2_0 + (1-xotd)* (1-pnumerator2_0))),        
                  den2=ifelse(plannedtrt==1,
                              (xotd*pdenominator2_1 + (1-xotd)* (1-pdenominator2_1)),
                              (xotd*pdenominator2_0 + (1-xotd)* (1-pdenominator2_0))),
                  num2= ifelse(time==0,1,num2),
                  den2= ifelse(time==0,1,den2)
    ) %>%
    dplyr:: group_by(id) %>%
    dplyr:: mutate(
      n2= cumprod(num2),
      d2= cumprod(den2)
    )%>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      stabwt2= n2/d2,
      wt2=1/d2
    )
  
  
  dd_wts_xo= dd[complete.cases(dd$xotd) & dd$xotd==0,]
  
  # Distribution of XO unstab wts
  mean_unstab_xo<- mean(dd_wts_xo$wt2)
  sd_unstab_xo<- sd(dd_wts_xo$wt2)
  min_unstab_xo<- min(dd_wts_xo$wt2)
  max_unstab_xo<- max(dd_wts_xo$wt2)
  
  # Distribution of XO stab wts
  mean_stab_xo<- mean(dd_wts_xo$stabwt2)
  sd_stab_xo<- sd(dd_wts_xo$stabwt2)
  min_stab_xo<- min(dd_wts_xo$stabwt2)
  max_stab_xo<- max(dd_wts_xo$stabwt2)
  
  
  
  ## Step 4: Remove observations at censoring for XO and ltfu and fit the outcomes model---------
  
  dd_outcome<- dd[dd$xotd==0 & dd$censoringtd_new==0,]
  
  
  # Cox proportional hazard model
  xo_lofu_cox <- coxph(Surv(time,time+1,deathtd) ~ plannedtrt,
                       data = dd_outcome, weights= stabwt1*stabwt2, ties = ("breslow"))
  ### HR
  HR<- exp(coef(xo_lofu_cox))[1]
  
  
  # Pooled logistic regression
  xo_lofu_glm<- glm(deathtd ~ time + time2 + factor(plannedtrt),
                    data= dd_outcome, weights= stabwt1*stabwt2, family =quasibinomial())
  ### dHR
  dHR<- exp(coef(xo_lofu_glm))[[4]]
  ### cHR and RR
  dd_outcome$prob<- 1-predict(xo_lofu_glm, newdata=dd_outcome, type="response")
  dd_outcome_1<- dd_outcome%>%
    dplyr::arrange(id,time)%>%
    dplyr::group_by(id)%>%
    dplyr::mutate(
      surv=cumprod(prob))
  dd_outcome_1<- dd_outcome_1[,c('surv','plannedtrt','time')]
  results<- dd_outcome_1%>%
    dplyr::group_by(time, plannedtrt)%>%
    dplyr::summarize(mean_survival=mean(surv))
  results<- results%>%
    dplyr::ungroup()%>%
    dplyr:: mutate(
      time=time+1
    )
  results<- dplyr::bind_rows(c(time=0,plannedtrt=0,mean_survival=1),c(time=0,plannedtrt=1,mean_survival=1),results )
  results$plannedtrtf<- factor(results$plannedtrt,labels=c("Crizotinib","Brigatinib"))
  wide_results<- reshape2::dcast(results, time ~ plannedtrtf,value.var = 'mean_survival')
  wide_results<- wide_results%>%
    mutate(
      risk0=(1-Crizotinib),
      risk1=(1-Brigatinib),
      RD= (1-Brigatinib)-(1-Crizotinib),
      cHR= log(Brigatinib)/log(Crizotinib),
      RR=(1-Brigatinib)/(1-Crizotinib),
      s0=Crizotinib,
      s1=Brigatinib)

  
  # KM Estimator
  km.fit.xo.lofu = survfit(Surv(time,time+1,deathtd) ~ plannedtrt, data = dd_outcome, weights = stabwt1*stabwt2)
  ### cHR
  S0<- summary(km.fit.xo.lofu, times=k, extend=T)$surv[1]
  S1<- summary(km.fit.xo.lofu, times=k, extend=T)$surv[2]
  cHR_km<- log(S1)/log(S0)
  ### RR
  R0_km<- 1-summary(km.fit.xo.lofu, times = k)$surv[1]
  R1_km<- 1-summary(km.fit.xo.lofu, times = k)$surv[2]
  RD_km<- R1_km-R0_km
  RR_km<- R1_km/R0_km
  
  ### Surv. prob in control
  km_p0_12<- summary(km.fit.xo.lofu, times=12, extend=T)$surv[1]
  km_p0_24<- summary(km.fit.xo.lofu, times=24, extend=T)$surv[1]
  km_p0_36<- summary(km.fit.xo.lofu, times=36, extend=T)$surv[1]
  # 
  ### Surv. prob in exp
  km_p1_12<- summary(km.fit.xo.lofu, times=12, extend=T)$surv[2]
  km_p1_24<- summary(km.fit.xo.lofu, times=24, extend=T)$surv[2]
  km_p1_36<- summary(km.fit.xo.lofu, times=36, extend=T)$surv[2]
  # 
  ## Return results-------------------------------------------------------------
  
  #print(head(dd_outcome))
  #assign("data_ipcw_final",dd_outcome,envir = .GlobalEnv)
  #assign("results_final",results,envir = .GlobalEnv)

  return(c(HR,
           dHR, 
           wide_results$risk0[which(wide_results$time==k-1)],
           wide_results$risk1[which(wide_results$time==k-1)],
           wide_results$RD[which(wide_results$time==k-1)],
           wide_results$cHR[which(wide_results$time==k-1)],
           wide_results$RR[which(wide_results$time==k-1)],
           cHR_km,
           RR_km,
           wide_results$s0[which(wide_results$time==12)],
           wide_results$s0[which(wide_results$time==24)],
           wide_results$s0[which(wide_results$time==36)],
           wide_results$s1[which(wide_results$time==12)],
           wide_results$s1[which(wide_results$time==24)],
           wide_results$s1[which(wide_results$time==36)],
           mean_unstab,
           sd_unstab,
           min_unstab,
           max_unstab,
           mean_stab,
           sd_stab,
           min_stab,
           max_stab,
           mean_unstab_xo,
           sd_unstab_xo,
           min_unstab_xo,
           max_unstab_xo,
           mean_stab_xo,
           sd_stab_xo,
           min_stab_xo,
           max_stab_xo,
           km_p0_12,
           km_p0_24,
           km_p0_36,
           km_p1_12,
           km_p1_24,
           km_p1_36
           ))  
}

# Set replicate failures to NA
boot.robust<- function(data,indices){
  tryCatch({analyses.boot(data,indices)}, error= function(e){NA})
}

set.seed(7268)

log.results<- boot(data= df_censoring_xotd, statistic= boot.robust, R=20) 
log.results$t0

# CoxPH
HR<- round(log.results$t0[[1]],2)
HR
conf_interval_HR <- boot.ci(log.results, conf = 0.95, type = "perc", index = 1)
conf_interval_HR

# Pooled logistic regression
dHR<- round(log.results$t0[[2]],2)
dHR
conf_interval_dHR <- boot.ci(log.results, conf = 0.95, type = "perc", index = 2)
conf_interval_dHR
#
cHR<- round(log.results$t0[[6]],2)
cHR
conf_interval_cHR <- boot.ci(log.results, conf = 0.95, type = "perc", index = 6)
conf_interval_cHR
#
RR<- round(log.results$t0[[7]],2)
RR
conf_interval_RR <- boot.ci(log.results, conf = 0.95, type = "perc", index = 7)
conf_interval_RR

# KM Estimator
cHR_km<-round(log.results$t0[[8]],2)
cHR_km
conf_interval_cHR_km<- boot.ci(log.results, conf = 0.95, type = "perc", index = 8)
conf_interval_cHR_km
#
RR_km<-round(log.results$t0[[9]],2)
RR_km
conf_interval_RR_km<- boot.ci(log.results, conf = 0.95, type = "perc", index = 9)
conf_interval_RR_km

# Parametric survival prob in control
p0_12<- round(log.results$t0[[10]],2)
conf_interval_p0_12<-  boot.ci(log.results, conf = 0.95, type = "perc", index = 10)
p0_24<- round(log.results$t0[[11]],2)
conf_interval_p0_24<-  boot.ci(log.results, conf = 0.95, type = "perc", index = 11)
p0_36<-round(log.results$t0[[12]],2)
conf_interval_p0_36<-  boot.ci(log.results, conf = 0.95, type = "perc", index = 12)
#
# parametric survival prob in exp
p1_12<- round(log.results$t0[[13]],2)
conf_interval_p1_12<-boot.ci(log.results, conf = 0.95, type = "perc", index = 13) 
p1_24<- round(log.results$t0[[14]],2)
conf_interval_p1_24<- boot.ci(log.results, conf = 0.95, type = "perc", index = 14)
p1_36<-round(log.results$t0[[15]],2)
conf_interval_p1_36<- boot.ci(log.results, conf = 0.95, type = "perc", index = 15)


# Distribution of LTFU unstab wts
mean_unstab<- round(as.numeric(log.results$t0[[16]]),2)
sd_unstab<- round(as.numeric(log.results$t0[[17]]),2)
mean_sd_unstab<- paste0(mean_unstab,"(",sd_unstab,")")
min_unstab<- round(log.results$t0[[18]],2)
max_unstab<- round(log.results$t0[[19]],2)
min_max_unstab<- paste0(min_unstab,"(",max_unstab,")")

# Distribution of LTFU stab wts
mean_stab<- round(as.numeric(log.results$t0[[20]]),2)
sd_stab<- round(as.numeric(log.results$t0[[21]]),2)
mean_sd_stab<- paste0(mean_stab,"(",sd_stab,")")
min_stab<- round(as.numeric(log.results$t0[[22]]),2)
max_stab<- round(as.numeric(log.results$t0[[23]]),2)
min_max_stab<- paste0(min_stab,"(",max_stab,")")

# Distribution of XO unstab wts
mean_unstab_xo<- round(as.numeric(log.results$t0[[24]]),2)
sd_unstab_xo<- round(as.numeric(log.results$t0[[25]]),2)
mean_sd_unstab_xo<- paste0(mean_unstab_xo,"(",sd_unstab_xo,")")
min_unstab_xo<- round(as.numeric(log.results$t0[[26]]),2)
max_unstab_xo<- round(as.numeric(log.results$t0[[27]]),2)
min_max_unstab_xo<- paste0(min_unstab_xo,"(",max_unstab_xo,")")

# Distribution of XO stab wts
mean_stab_xo<- round(as.numeric(log.results$t0[[28]]),2)
sd_stab_xo<- round(as.numeric(log.results$t0[[29]]),2)
mean_sd_stab_xo<- paste0(mean_stab_xo,"(",sd_stab_xo,")")
min_stab_xo<- round(as.numeric(log.results$t0[[30]]),2)
max_stab_xo<- round(as.numeric(log.results$t0[[31]]),2)
min_max_stab_xo<- paste0(min_stab_xo,"(",max_stab_xo,")")
  

# KM survival prob in control
km_p0_12<- round(log.results$t0[[32]],2)
conf_interval_km_p0_12<-  boot.ci(log.results, conf = 0.95, type = "perc", index = 32)
km_p0_24<- round(log.results$t0[[33]],2)
conf_interval_km_p0_24<-  boot.ci(log.results, conf = 0.95, type = "perc", index = 33)
km_p0_36<-round(log.results$t0[[34]],2)
conf_interval_km_p0_36<-  boot.ci(log.results, conf = 0.95, type = "perc", index = 34)
#
# KM survival prob in exp
km_p1_12<- round(log.results$t0[[35]],2)
conf_interval_km_p1_12<-boot.ci(log.results, conf = 0.95, type = "perc", index = 35) 
km_p1_24<- round(log.results$t0[[36]],2)
conf_interval_km_p1_24<- boot.ci(log.results, conf = 0.95, type = "perc", index = 36)
km_p1_36<-round(log.results$t0[[37]],2)
conf_interval_km_p1_36<- boot.ci(log.results, conf = 0.95, type = "perc", index = 37)
#
#
#
# spec 4x4 results (20 bootstrap replicates)---------------------------------

spec_4x4_stab<- data.frame(
  Method=  c("HR","dHR","cHR","RR","cHR_km","RR_km"),
  Estimate= c(HR,dHR,cHR,RR,cHR_km,RR_km),
  LL= c(0.6595,0.6225,0.6269,0.7156,0.3989,0.5352),
  UL= c(1.4143,1.3893,1.3851,1.2974,1.2220,1.1648)
)
spec_4x4_stab[,sapply(spec_4x4_stab, is.numeric)]<- round(spec_4x4_stab[,sapply(spec_4x4_stab, is.numeric)],2)
write.csv(spec_4x4_stab,file="spec_4x4_stab.csv",row.names = T)
#
# spec 4x4 km prob (20 bootstrap replicates) -----------------------------------------------------
km_prob_ipcw_4x4 <- data.frame(
  Time=  c("12","24","36"),
  Estimate_control= c(km_p0_12,km_p0_24,km_p0_36),
  LL_control= c(0.9822,0.8990,0.5707),
  UL_control= c(1.0000,0.9767,0.8509),
  Estimate_exp= c(km_p1_12,km_p1_24,km_p1_36),
  LL_exp= c(0.9949,0.8981,0.6809),
  UL_exp= c(1.0000,0.9289,0.7173)
)
km_prob_ipcw_4x4[,sapply(km_prob_ipcw_4x4, is.numeric)]<- round(km_prob_ipcw_4x4[,sapply(km_prob_ipcw_4x4, is.numeric)],2)*100
write.csv(km_prob_ipcw_4x4,file="km_prob_ipcw_4x4.csv",row.names = T)

# KM Plot IPCW 4x4--------------------------------------------------------------

# Analysis for Figure 5 in the manuscript
# Note: run spec 4x4  or spec 4x7
# Change name or prob table and cHR and RR text accordingly
windowsFonts(x=windowsFont("Arial"))

# Spec 4 for LTFU
# Numerator control for LTFU
nmodel_lofu_0<- glm(censoringtd_new ~ time + time2,
                    data=df_ipcw_pp[df_ipcw_pp$plannedtrt==0,], family =binomial())

# Numerator exp for LTFU
nmodel_lofu_1<- glm(censoringtd_new ~ time + time2, 
                    data=df_ipcw_pp[df_ipcw_pp$plannedtrt==1,], family =binomial())

# Predict
df_ipcw_pp$pnumerator1_0<- predict(nmodel_lofu_0, newdata=df_ipcw_pp, type='response')
df_ipcw_pp$pnumerator1_1<- predict(nmodel_lofu_1, newdata=df_ipcw_pp, type='response')

# Denominator control for LTFU
dmodel_lofu_0<- glm(
  censoringtd_new ~ time + SMKHISGR1 + STRATRAN + ecogtdGR1
  + AGE + MICNSFL + tltd_wins + progtd + trttd + LIASE + DIAGINITGR3
  + time2
  + AGE2,
  data=df_ipcw_pp[df_ipcw_pp$plannedtrt==0,],
  family =binomial())

# Denominator exp for LTFU
dmodel_lofu_1<- glm(
  censoringtd_new ~ time + AGE + SEX + LIASE + STRATRAN + ECOGGR1 + ecogtdGR1 
  + tltd_wins + icprogtd + progtd + RACEGR1 + SMKHISGR1 + MICNSFL 
  + DIAGINITGR3 + PRRADYN + time2 + AGE2,
  data=df_ipcw_pp[df_ipcw_pp$plannedtrt==1,],
  family =binomial())

# Predict
df_ipcw_pp$pdenominator1_0<- predict(dmodel_lofu_0, newdata=df_ipcw_pp, type='response')
df_ipcw_pp$pdenominator1_1<- predict(dmodel_lofu_1, newdata=df_ipcw_pp, type='response')


# Weight Calculation

# Order the data by "id" and "time"
df_ipcw_pp<- df_ipcw_pp[order(df_ipcw_pp$id,df_ipcw_pp$time),]

# Calculate num1 and den1 
df_ipcw_pp<- df_ipcw_pp%>%dplyr::mutate(
  num1= ifelse(plannedtrt==1,
               (censoringtd_new*pnumerator1_1 + (1-censoringtd_new)*(1-pnumerator1_1)),
               (censoringtd_new*pnumerator1_0 + (1-censoringtd_new)* (1-pnumerator1_0))),
  den1=ifelse(plannedtrt==1,
              (censoringtd_new*pdenominator1_1 + (1-censoringtd_new)* (1-pdenominator1_1)),
              (censoringtd_new*pdenominator1_0 + (1-censoringtd_new)* (1-pdenominator1_0))),
  num1= ifelse(time==0,1,num1),
  den1= ifelse(time==0,1,den1)
) %>%
  dplyr::group_by(id) %>%
  dplyr::mutate(
    n1= cumprod(num1),
    d1= cumprod(den1)
  )%>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    stabwt1= n1/d1, # stabilized weights
    wt1=1/d1        # unstabilized weights 
  )

# spec 4
# 
### Numerator control for XO
nmodel_xotd_0<- glm( xotd ~ time, data=df_ipcw_pp[complete.cases(df_ipcw_pp$xotd) & df_ipcw_pp$plannedtrt==0,], family =binomial()) # Predict P for control arm

# Predict
df_ipcw_pp$pnumerator2_0<- predict(nmodel_xotd_0, newdata=df_ipcw_pp, type='response') #Set Numerator for Experimental =1 by setting P to 0
df_ipcw_pp$pnumerator2_1<- 0


### Denominator control for XO
dmodel_xotd_0<- glm(
  xotd ~ time
  + SMKHISGR1 + ecogtdGR1 + AGE + tltd_wins + progtime
  + STRATRAN +  LIASE + DIAGINITGR3 + MICNSFL
  + I(progtime^2),
  data=df_ipcw_pp[complete.cases(df_ipcw_pp$xotd) & df_ipcw_pp$plannedtrt==0 & df_ipcw_pp$progtd==1,],
  family =binomial())

# Predict
df_ipcw_pp$pdenominator2_0<- predict(dmodel_xotd_0, newdata=df_ipcw_pp, type='response') # Set Denominator for Control non-progressors to 1 by setting P to 0
df_ipcw_pp$pdenominator2_0[df_ipcw_pp$plannedtrt==0 & df_ipcw_pp$progtd==0]=0 # Set Denominator for Experimental to 1 by setting P to 0
df_ipcw_pp$pdenominator2_1<- 0

# spec 7

# ### Numerator control for XO
# nmodel_xotd_0<- glm( xotd ~ time, data=df_ipcw_pp[complete.cases(df_ipcw_pp$xotd) & df_ipcw_pp$plannedtrt==0,],
#                      family =binomial()) # Predict P for control arm
# 
# # Predict
# df_ipcw_pp$pnumerator2_0<- predict(nmodel_xotd_0, newdata=df_ipcw_pp, type='response') # Set Numerator for Experimental =1 by setting P to 0
# df_ipcw_pp$pnumerator2_1<- 0
# 
# 
# ### Denominator control for XO
# dmodel_xotd_0<- glm(
#   xotd ~ time
#   + SMKHISGR1 + ecogtdGR2 + AGE + tltd_wins
#   + STRATRAN +  LIASE + DIAGINITGR3 + MICNSFL,
#   data=df_ipcw_pp[complete.cases(df_ipcw_pp$xotd) & df_ipcw_pp$plannedtrt==0 & df_ipcw_pp$progtd==1,],
#   family =binomial())
# 
# # Predict
# df_ipcw_pp$pdenominator2_0<- predict(dmodel_xotd_0, newdata=df_ipcw_pp, type='response') # Set Denominator for Control non-progressors to 1 by setting P to 0
# df_ipcw_pp$pdenominator2_0[df_ipcw_pp$plannedtrt==0 & df_ipcw_pp$progtd==0]=0 # Set Denominator for Experimental to 1 by setting P to 0
# df_ipcw_pp$pdenominator2_1<- 0
# 

# Weight Calculation

df_ipcw_pp<- df_ipcw_pp[order(df_ipcw_pp$id,df_ipcw_pp$time),] # Contribution from xotd(t)=0 is 1-P # Contribution from xotd(t)=1 is P
df_ipcw_pp<- df_ipcw_pp%>%
  dplyr::mutate(num2= ifelse(plannedtrt==1,
                             (xotd*pnumerator2_1 + (1-xotd)*(1-pnumerator2_1)),
                             (xotd*pnumerator2_0 + (1-xotd)* (1-pnumerator2_0))),        
                den2=ifelse(plannedtrt==1,
                            (xotd*pdenominator2_1 + (1-xotd)* (1-pdenominator2_1)),
                            (xotd*pdenominator2_0 + (1-xotd)* (1-pdenominator2_0))),
                num2= ifelse(time==0,1,num2),
                den2= ifelse(time==0,1,den2)
  ) %>%
  dplyr:: group_by(id) %>%
  dplyr:: mutate(
    n2= cumprod(num2),
    d2= cumprod(den2)
  )%>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    stabwt2= n2/d2,
    wt2=1/d2
  )

dd_outcome<- df_ipcw_pp[df_ipcw_pp$xotd==0 & df_ipcw_pp$censoringtd_new==0,]

kmfit_final = survfit(Surv(time,time+1,deathtd) ~ plannedtrt, data = dd_outcome,
                 weights = stabwt1*stabwt2)
survplot_ipcw<- ggsurvplot(kmfit_final,
                           data=dd_outcome,
                           #title="Overall Survival: IPCW adjusted analysis for switching and \n Loss to follow-up",
                           ylab="Overall survival (%)",
                           xlab="Time after randomisation (months)",
                           size=0.6,
                           conf.int = F,
                           censor=T,
                           #censor.shape="|",
                           censor.size=3,
                           ncensor.plot=F,
                           legend.title="",
                           legend.labs= c("Crizotinib","Brigatinib"),
                           xlim=c(0,36),
                           ylim= c(0,100),
                           break.time.by=6,
                           ggtheme = theme_classic2(base_size = 14, base_family = "x"),
                           font.family="x",    # size of base 
                           risk.table = F,     #"nrisk_cumevents",
                           fun="pct",          # plot cumulative incidence,
                           surve.scale='percent',
                           cumevents = F,
                           risk.table.height=0.13,
                           risk.table.y.text.col=T,
                           risk.table.y.text=F,
                           fontsize=4,          # data size in table
                           font.title=9,        # size of table title
                           legend=c(0.15,0.65), # position of legend
                           palette = c("#00539CFF","#F95700FF"),
                           tables.theme = clean_theme())

# plot.margin = c(top, right, bottom, left)
survplot_ipcw$plot<- survplot_ipcw$plot +
  theme(plot.margin = unit(c(10,20,10,30),"points"))

#plot.margin = c(top, right, bottom, left)
# survplot_ipcw$table<- survplot_ipcw$table +
#   theme(
#     plot.title = element_text(size=12,color="black", hjust=-0.049),
#     plot.margin = unit(c(2,20,0,30),"points"))

# add text to plot
survplot_ipcw$plot<- survplot_ipcw$plot +
  ggplot2::annotate(
    "text",
    x=Inf, y=Inf,
    vjust=6, hjust=2.3,
    label="cHR 0.90 (95% CI 0.63,1.39)\n RR 0.92 (95% CI 0.72,1.30)",
    size=4)

survplot_ipcw # Save as pdf: 6x8
#
#
#
# Parametric plot IPCW 4x4 -----------------------------------------------------

# Analysis for Figure 7a in the manuscript

windowsFonts(x=windowsFont("Arial"))


# Numerator control for LTFU
nmodel_lofu_0<- glm(censoringtd_new ~ time + time2,
                    data=df_ipcw_pp[df_ipcw_pp$plannedtrt==0,], family =binomial())

# Numerator exp for LTFU
nmodel_lofu_1<- glm(censoringtd_new ~ time + time2, 
                    data=df_ipcw_pp[df_ipcw_pp$plannedtrt==1,], family =binomial())

# Predict
df_ipcw_pp$pnumerator1_0<- predict(nmodel_lofu_0, newdata=df_ipcw_pp, type='response')
df_ipcw_pp$pnumerator1_1<- predict(nmodel_lofu_1, newdata=df_ipcw_pp, type='response')

# Denominator control for LTFU
dmodel_lofu_0<- glm(
  censoringtd_new ~ time + SMKHISGR1 + STRATRAN + ecogtdGR1
  + AGE + MICNSFL + tltd_wins + progtd + trttd + LIASE + DIAGINITGR3
  + time2
  + AGE2,
  data=df_ipcw_pp[df_ipcw_pp$plannedtrt==0,],
  family =binomial())

# Denominator exp for LTFU
dmodel_lofu_1<- glm(
  censoringtd_new ~ time + AGE + SEX + LIASE + STRATRAN + ECOGGR1 + ecogtdGR1 
  + tltd_wins + icprogtd + progtd + RACEGR1 + SMKHISGR1 + MICNSFL 
  + DIAGINITGR3 + PRRADYN + time2 + AGE2,
  data=df_ipcw_pp[df_ipcw_pp$plannedtrt==1,],
  family =binomial())

# Predict
df_ipcw_pp$pdenominator1_0<- predict(dmodel_lofu_0, newdata=df_ipcw_pp, type='response')
df_ipcw_pp$pdenominator1_1<- predict(dmodel_lofu_1, newdata=df_ipcw_pp, type='response')


# Weight Calculation

# Order the data by "id" and "time"
df_ipcw_pp<- df_ipcw_pp[order(df_ipcw_pp$id,df_ipcw_pp$time),]

# Calculate num1 and den1 
df_ipcw_pp<- df_ipcw_pp%>%dplyr::mutate(
  num1= ifelse(plannedtrt==1,
               (censoringtd_new*pnumerator1_1 + (1-censoringtd_new)*(1-pnumerator1_1)),
               (censoringtd_new*pnumerator1_0 + (1-censoringtd_new)* (1-pnumerator1_0))),
  den1=ifelse(plannedtrt==1,
              (censoringtd_new*pdenominator1_1 + (1-censoringtd_new)* (1-pdenominator1_1)),
              (censoringtd_new*pdenominator1_0 + (1-censoringtd_new)* (1-pdenominator1_0))),
  num1= ifelse(time==0,1,num1),
  den1= ifelse(time==0,1,den1)
) %>%
  dplyr::group_by(id) %>%
  dplyr::mutate(
    n1= cumprod(num1),
    d1= cumprod(den1)
  )%>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    stabwt1= n1/d1, # stabilized weights
    wt1=1/d1        # unstabilized weights 
  )


### Numerator control for XO
nmodel_xotd_0<- glm( xotd ~ time, data=df_ipcw_pp[complete.cases(df_ipcw_pp$xotd) & df_ipcw_pp$plannedtrt==0,], family =binomial()) # Predict P for control arm

# Predict
df_ipcw_pp$pnumerator2_0<- predict(nmodel_xotd_0, newdata=df_ipcw_pp, type='response') #Set Numerator for Experimental =1 by setting P to 0
df_ipcw_pp$pnumerator2_1<- 0


### Denominator control for XO
dmodel_xotd_0<- glm(
  xotd ~ time
  + SMKHISGR1 + ecogtdGR1 + AGE + tltd_wins + progtime
  + STRATRAN +  LIASE + DIAGINITGR3 + MICNSFL
  + I(progtime^2),
  data=df_ipcw_pp[complete.cases(df_ipcw_pp$xotd) & df_ipcw_pp$plannedtrt==0 & df_ipcw_pp$progtd==1,],
  family =binomial())

# Predict
df_ipcw_pp$pdenominator2_0<- predict(dmodel_xotd_0, newdata=df_ipcw_pp, type='response') # Set Denominator for Control non-progressors to 1 by setting P to 0
df_ipcw_pp$pdenominator2_0[df_ipcw_pp$plannedtrt==0 & df_ipcw_pp$progtd==0]=0 # Set Denominator for Experimental to 1 by setting P to 0
df_ipcw_pp$pdenominator2_1<- 0



# Weight Calculation

df_ipcw_pp<- df_ipcw_pp[order(df_ipcw_pp$id,df_ipcw_pp$time),] # Contribution from xotd(t)=0 is 1-P # Contribution from xotd(t)=1 is P
df_ipcw_pp<- df_ipcw_pp%>%
  dplyr::mutate(num2= ifelse(plannedtrt==1,
                             (xotd*pnumerator2_1 + (1-xotd)*(1-pnumerator2_1)),
                             (xotd*pnumerator2_0 + (1-xotd)* (1-pnumerator2_0))),        
                den2=ifelse(plannedtrt==1,
                            (xotd*pdenominator2_1 + (1-xotd)* (1-pdenominator2_1)),
                            (xotd*pdenominator2_0 + (1-xotd)* (1-pdenominator2_0))),
                num2= ifelse(time==0,1,num2),
                den2= ifelse(time==0,1,den2)
  ) %>%
  dplyr:: group_by(id) %>%
  dplyr:: mutate(
    n2= cumprod(num2),
    d2= cumprod(den2)
  )%>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    stabwt2= n2/d2,
    wt2=1/d2
  )

dd_outcome<- df_ipcw_pp[df_ipcw_pp$xotd==0 & df_ipcw_pp$censoringtd_new==0,]
                        
# Pooled logistic regression
xo_lofu_glm<- glm(deathtd ~ time + time2 + factor(plannedtrt),
                  data= dd_outcome, weights= stabwt1*stabwt2, family =quasibinomial())

dd_outcome$prob<- 1-predict(xo_lofu_glm, newdata=dd_outcome, type="response")
dd_outcome_1<- dd_outcome%>%
  dplyr::arrange(id,time)%>%
  dplyr::group_by(id)%>%
  dplyr::mutate(
    surv=cumprod(prob))
dd_outcome_1<- dd_outcome_1[,c('surv','plannedtrt','time')]
results<- dd_outcome_1%>%
  dplyr::group_by(time, plannedtrt)%>%
  dplyr::summarize(mean_survival=mean(surv))
results<- results%>%
  dplyr::ungroup()%>%
  dplyr:: mutate(
    time=time+1
  )
results<- dplyr::bind_rows(c(time=0,plannedtrt=0,mean_survival=1),c(time=0,plannedtrt=1,mean_survival=1),results )
results$plannedtrtf<- factor(results$plannedtrt,labels=c("Crizotinib","Brigatinib"))

results_plot<- results
results_plot$mean_survival<- results_plot$mean_survival*100
results_plot<- na.omit(results_plot)

ggplot2<- ggplot(
  results_plot,aes(x=time,y=mean_survival))+
  geom_line(aes(colour=plannedtrtf))+
  geom_point(aes(colour=plannedtrtf))+
  xlab("Months")+
  scale_x_continuous(limits=c(0,36),breaks=seq(0,36,6))+
  ylab("Overall Survival (%)")+
  scale_y_continuous(limits=c(0,100),breaks=seq(0,100, by=25))+
  #ggtitle("IPCW: Standardized Survival Curve Weighted for Baseline and\n Time Varying Confounders")+
  labs(colour="Treatment")+
  theme_bw()+
  theme(legend.position = c(0.2,0.6),
        panel.border = element_blank(),
        legend.background = element_rect(fill="transparent"),
        #legend.key= element_rect(fill="transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color="black"))+
  scale_color_manual(values=c("#00539CFF","#F95700FF"))
ggplot2

ggplot2+theme(axis.text = element_text(size=15),
              axis.title = element_text(size=15),
              #title=element_text(size=10),
              #legend.title=element_text(size=14),
              legend.text=element_text(size=10))

# add text to plot
ggplot2 + annotate("text", x=20,y=70, label="cHR 0.90 (95% CI 0.63,1.39)\n RR 0.92 (95% CI 0.72,1.30)", size=unit(4,"pt"))


## END------
