
# File name: gformula_analysis_synthetic.R
# Author: Amani Al Tawil
# Date: 25th January 2024
# Data set: df_ipcw_pp (synthetic data in person-time format) 
# Method: Parametric g-formula analyses to adjust for LTFU & switching 
# References: McGrath 2020 
#
#
#
# Description:
#
# This R script performs parametric g-Formula analyses on the generated synthetic data
# to account for loss to follow-up (LTFU) and treatment switching (XO). 
# The analysis is based on the methodology described by McGrath et al 2020
# Refer to Electronic Supplementary Material (ESM) 1 for the data dictionary

# We performed five different analyses with different model specifications for the time-varying covariates (TVC) and outcome. 
## In specification 1, for each fitted model, we considered including all predictors with a P-value < 0.2. 
## In specifications 2 and 3, for each fitted model, we considered including all predictors with a P-value < 0.1 and 0.05, respectively. 
## In Specification 4, we repeated specification 1, but introduced three levels for the time-varying ECOG variable instead of two levels. 
## In specification 5, we repeated specification 1, but  removed time to progression "progtime" variable from the outcome model.
## Specifications 4 and 5 were included to enable meaningful comparisons with other specifications from the IPCW approach
## Refer to ESM 4 for the list of Covariates included in each specification of the parametric g-formula analysis
## Since we only provide this code applied to the synthetic data,
## we do not include details regarding our model selection process since it was applied on the actual data
## and is therefore irrelevant here

# Usage:
#
# 1. Set the working directory to the location of your data set and this script.
# 2. Ensure that the required libraries are installed and loaded.
# 3. Customize the parameters and variables as needed for your specific analysis.
# 4. Run the script to perform the parametric g-Formula analysis.


# Notes: 
#
# For details about the functionality of the package reference to McGrath 2020 and
# corresponding supplementary material
# To access the code that is used to generate the synthetic data, refer to 'synthetic_data.R'
# For guidance on how to structure your data to be applied to the gfoRmula package refer to McGrath 2020
#
#
#
# Set Working directory---------------------------------------------------------

# Replace "YOUR_DIRECTORY_PATH" with the actual directory path.
setwd("YOUR_DIRECTORY_PATH") 
#setwd("C:\\Users\\altawil\\Desktop\\myproject_github")
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
#
#
#
# Install required Packages-----------------------------------------------------

# This will automatically install them into your default library. 

if(!require("dplyr"))install.packages("dplyr")
library(dplyr)
if(!require("survival"))install.packages("survival")
library(survival)
if(!require("survminer"))install.packages("survminer")
library(survminer)
if(!require("ggplot2"))install.packages("ggplot2")
library(ggplot2)
if(!require("openxlsx"))install.packages("openxlsx")
library(openxlsx)
if(!require("data.table"))install.packages("data.table")
library(data.table)
if(!require("forestplot"))install.packages("forestplot")
library(forestplot)


# gfoRmula Package Version:
#
# Our analysis uses the development version of the gfoRmula package, specifically
news(package= "gfoRmula") # Changes in version 1.0.4."
#
# To Install the Development Version:
# 1. Ensure you have the 'devtools' package installed by running: install.packages("devtools")
# 2. Then, you can install the development version of gfoRmula from GitHub using:
devtools::install_github("CausalInference/gfoRmula")
#
# For CRAN Version (not used for this analysis):
# Uncomment the following line if you wish to load the CRAN version of gfoRmula package.
#if (!require("gfoRmula")) install.packages("gfoRmula") 
# Please note that this analysis specifically relies on the GitHub version due to specific requirements.
#
#
#
# Parametric g-formula analyses:------------------------------------------------
#
#
#
# Analysis:
#
# Commenting the code:
# In the following section, lines 121-413 entitled "Specification 1: Apply gformula in Control and Exp arm (full model)"
# We provide comments to explain how the code works in our specific context
# To avoid redundancy, we only comment specification 1 analysis 
# The comments are applicable to subsequent specifications (2 to 5)
# Please review and adjust the code to fit your specific data and research context.
#
#
#
#
# Specification 1: Apply gformula in Control and Exp arm (full model) ----------
# 
#
#
# Specification 1: Apply gformula in control


# Data set:
# Load the data set "df_ipcw_pp" for this analysis 
# Subset the data to include only control group as we perform this analysis separately for each arm
# Convert the final.gformula_control to a data.table 
df_gformula<- df_ipcw_pp 
final.gformula_control<- subset(df_gformula, df_gformula$plannedtrt==0) 
final.gformula_control<- as.data.table(final.gformula_control)


# Custom history and function explanation:
# This custom history function is designed for the "progtime" variable, representing time to progression.
# The "progtime" value depends on changes in the TVC progression "progtd"
# It cannot be treated as a baseline or modeled as TVC.
# This function calculates "progtime" based on the state of "progtd" at different time points.
# At t = 0, it initializes 'progtime' for individuals: progtime=0 if progtd=1  or NA if progtd=0
# For subsequent time points (t > 0), it updates "progtime" based on the state of progression.
# Note: Customize this function according to your specific needs.

time_since_switch<- function(pool, histvars, time_name, t,id_name){
  if(t== 0){
    pool[progtd==1 & get(time_name)==t,':=' (progtime=0, now.progtd=1)]
    pool[progtd==0 & get(time_name)==t,':=' (progtime=NA, now.progtd=NA)]
  } else{
    pool[progtd==1 & lag1_progtd==0 & get(time_name)==t,':=' (progtime=t, now.progtd=1)]
    pool[progtd==1 & lag1_progtd==1 & get(time_name)==t,':=' (now.progtd=0)]
  }
}

# Set up parameters and variables for the gformula analysis:
obs_data = final.gformula_control
class(obs_data)
outcome_name = 'deathtd'
outcome_type = 'survival'
censor_name= 'censoringtd_new'
id='id'
time_points= 47 
time_name= 'time'
covnames = c('progtd','ecogtdGR1','icprogtd','tltd_wins','trttd')    # List your TVCs here...
covtypes = c('absorbing','binary','absorbing','normal','absorbing')  # Define the type of TVCs

# History and respective variables of interest:
histories = c(lagged,time_since_switch) # Define the history of TVC of interest (other options are cumulative lag, cumulative average lag)                        
histvars = list(c('progtd','ecogtdGR1','icprogtd','tltd_wins','trttd'),c('progtd')) # Define the TVC for which you are interested to generate histories

# Model specifications:
# Define the different models for the TVCs and treatment
# Note: Use only the lag value of the TVC if the variable has not been previously modeled 
# For e.g. you can only use icprogtd as predictor in the tltd_wins after it has been modeled
# You cannot use icprogtd as predictor in the progtd and ecogtdGR1 model as it has not been defined 
covparams= list(covmodels= c(
  progtd ~ time + I(time^2) + RACEGR1 + ECOGGR1 + SMKHISGR1 + lag1_icprogtd + lag1_tltd_wins, 
  
  ecogtdGR1 ~ time + AGE + I(AGE^2) + SEX + RACEGR1 + ECOGGR1 + SMKHISGR1+
    MICNSFL + PRRADYN + STRATRAN  + lag1_tltd_wins  + lag1_icprogtd + lag1_progtd + progtd + lag1_trttd,
  
  icprogtd ~ time + AGE + RACEGR1 + PRRADYN + STRATRAN + lag1_progtd + progtd + lag1_tltd_wins + lag1_trttd,
  
  tltd_wins ~  time + I(time^2) + AGE  + SEX + RACEGR1 + ECOGGR1 + MICNSFL + PRRADYN + STRATRAN +
    lag1_ecogtdGR1 +lag1_icprogtd +lag1_progtd +ecogtdGR1 + progtd + icprogtd + lag1_trttd,
  # Here we can comment the trttd model, since our regimen of interest is static 
  trttd ~  time + AGE + I(AGE^2) + RACEGR1 + STRATRAN  + lag1_icprogtd + icprogtd +  progtime + I(progtime^2))) 

# Define the outcome model:
ymodel <- deathtd ~ time + I(time^2) + AGE + I(AGE^2) + SMKHISGR1  + STRATRAN  + lag1_ecogtdGR1 +
  ecogtdGR1  + lag1_trttd + trttd + progtime

# Define the censoring model (Optional, could be commented out):
# It could be commented out since it is also treated as a static intervention
# We can fit a censoring model, if we are interested in computing IPCW estimates
censor_model<- censoringtd_new ~ time + I(time^2) + AGE +  I(AGE^2)+ SEX + RACEGR1+
  ECOGGR1 + SMKHISGR1  + PRRADYN + STRATRAN  + lag1_icprogtd + icprogtd +
  lag1_tltd_wins + tltd_wins + lag1_progtd + progtd + ecogtdGR1 + lag1_ecogtdGR1 +
  lag1_trttd + trttd 

# Treatment Variable or Intervention:
intvars= list('trttd')
interventions <- list(list(c(static, rep(0,time_points)))) # Type of intervention: Here it is static which sets treatment "trttd" to 0 at all time points
int_descript <- c('Never treat') 

# Nb of simulations and bootstrap samples:
#
# For computing CIs around the effect estimates,
# it is required to set the number of bootstrap samples using
# nsamples= x where x>0 
# For the sake of demonstration and to limit the computational time,
# We used nsamples= 20 here.
# In your actual data analysis,consider using a larger number of samples
# e.g. nsamples=500 or 1000
  
nsimul <- 1100 # Set to a larger value than the number of observation
nsamples<- 20  # Set nsamples > 0 to estimate percentile based bootstrap CI 

# Parallel Processing Configuration:
ncores <- parallel::detectCores()-1

# Run the gformula
gform_1 <- gformula(obs_data = obs_data,
                    outcome_type = outcome_type,
                    id=id,
                    time_points = time_points,
                    time_name = time_name,
                    covnames = covnames,
                    outcome_name = outcome_name,
                    covtypes = covtypes, covparams = covparams,
                    ymodel = ymodel, intvars = intvars,
                    interventions = interventions,
                    int_descript = int_descript,histories = histories,
                    histvars = histvars,
                    nsimul = nsimul,
                    seed = 1234,
                    nsamples= nsamples,
                    parallel=TRUE,
                    ncores=ncores,
                    ci_method = 'percentile',
                    boot_diag = TRUE,
                    model_fits = TRUE,
                    censor_name=censor_name,   # could be commented out
                    censor_model=censor_model, # could be commented out
                    ipw_cutoff_quantile = 0.95,
                    # Define baseline covariates:
                    basecovs = c('AGE','SEX','RACEGR1','ECOGGR1','SMKHISGR1','MICNSFL','PRRADYN','STRATRAN'),
                    # Incorporate deterministic information in the form of restrictions to help mitigate bias from model misspecification
                    restrictions = list(c('trttd','progtd==1',simple_restriction,0)), # In our case time-varying treatment "trttd"= 0 when progression status "progtd" =0  in the control arm
                    # Incorporate deterministic information in the form of outcome restrictions:
                    yrestrictions = list(c('progtd==1',0)),# In our case "deathtd"=0 when "progression"progtd" != 1 
                    # Set sim_data_b = F with nsamples>0 and sim_data_b = T with nsamples =0 to inspect the simulated data
                    sim_data_b = F)

# Display the Results:
print(gform_1)
#
#
#
# Inspect simulated data
# Rerun gform_1 with nsamples=0 and sim_data=T
sim_nc<- gform_1$sim_data$'Natural course'
table(sim_nc$progtd, sim_nc$trttd)
table(sim_nc$progtd, sim_nc$Y)
table(sim_nc$progtd, sim_nc$now.progtd,useNA = "always")
sim_nt<- gform_1$sim_data$`Never treat`
table(sim_nt$progtd, sim_nt$trttd)
table(sim_nt$progtd, sim_nt$Y)
#
#
#
#

# Specification 1: Apply gformula in Exp arm

# Subset the data to include only experimental group as we perform this analysis separately for each arm
df_gformula<- df_ipcw_pp
final.gformula_exp<- subset(df_gformula, df_gformula$plannedtrt==1)
final.gformula_exp<- as.data.table(final.gformula_exp)
#
# Custom history
time_since_switch<- function(pool, histvars, time_name, t,id_name){
  if(t== 0){
    pool[progtd==1 & get(time_name)==t,':=' (progtime=0, now.progtd=1)]
    pool[progtd==0 & get(time_name)==t,':=' (progtime=NA, now.progtd=NA)]
  } else{
    pool[progtd==1 & lag1_progtd==0 & get(time_name)==t,':=' (progtime=t, now.progtd=1)]
    pool[progtd==1 & lag1_progtd==1 & get(time_name)==t,':=' (now.progtd=0)]
  }
}
#
obs_data = final.gformula_exp
class(obs_data)
outcome_name = 'deathtd'
outcome_type = 'survival'
censor_name= 'censoringtd_new'
id='id'
time_points= 47
time_name= 'time'
covnames = c('progtd','ecogtdGR1','icprogtd','tltd_wins','trttd')
covtypes = c('absorbing','binary','absorbing','normal','binary')
histories = c(lagged,time_since_switch)
histvars = list(c('progtd','ecogtdGR1','icprogtd','tltd_wins','trttd'),c('progtd'))
covparams= list(covmodels= c(progtd ~ time + SEX + RACEGR1 + ECOGGR1 + SMKHISGR1 + STRATRAN +
                               lag1_icprogtd + lag1_tltd_wins,
                             
                             ecogtdGR1 ~ time + I(time^2) + AGE + I(AGE^2) +SEX + RACEGR1 + ECOGGR1 +SMKHISGR1+
                               MICNSFL + PRRADYN + STRATRAN  +
                               lag1_tltd_wins + lag1_icprogtd , 

                             icprogtd ~ time + I(time^2) + ECOGGR1 + STRATRAN  + progtd +
                               lag1_progtd + lag1_tltd_wins ,
                             
                             tltd_wins ~ time + I(time^2) + AGE +  I(AGE^2) + SEX 
                             + RACEGR1 + ECOGGR1 + STRATRAN + ecogtdGR1 + progtd + lag1_ecogtdGR1 +
                               lag1_progtd + lag1_icprogtd + icprogtd,
                            
                             trttd ~ 1))  # Since no switching in experimental arm, trttd is always 1

ymodel <- deathtd ~  time  + SEX + MICNSFL + STRATRAN + ecogtdGR1 +
  icprogtd  +lag1_ecogtdGR1 + lag1_icprogtd + progtime

censor_model<- censoringtd_new ~ time + I(time^2)+ ECOGGR1 + SMKHISGR1 +
  ecogtdGR1 + icprogtd  + progtd + lag1_ecogtdGR1 + lag1_icprogtd + lag1_progtd

intvars= list('trttd')
interventions <- list(list(c(static, rep(1,time_points))))
int_descript <- c('Always treat')
nsimul <- 1100
nsamples<- 20 # 500
ncores <- parallel::detectCores()-1
#
gform_2 <- gformula(obs_data = obs_data,
                    outcome_type = outcome_type,
                    id=id,
                    time_points = time_points,
                    time_name = time_name,
                    covnames = covnames,
                    outcome_name = outcome_name,
                    covtypes = covtypes, 
                    covparams = covparams,
                    ymodel = ymodel,
                    intvars = intvars,
                    interventions = interventions,
                    int_descript = int_descript,
                    histories = histories,
                    histvars = histvars,
                    basecovs = c('AGE','SEX','RACEGR1','ECOGGR1','SMKHISGR1','MICNSFL','PRRADYN','STRATRAN'),                    
                    yrestrictions = list(c('progtd==1',0)),
                    nsimul = nsimul,
                    seed = 1234,
                    nsamples= nsamples,
                    parallel=TRUE,
                    ncores=ncores,
                    ci_method = 'percentile',
                    censor_name=censor_name,
                    censor_model=censor_model,
                    ipw_cutoff_quantile = 0.95,#
                    boot_diag = TRUE,
                    model_fits = TRUE)
print(gform_2)
#
#
#
#### Specification 1: RR and cHR------------------------------------------------

# Estimate Relative Risk (RR) as the ratio of g-risks under the two treatment strategies
# RR= g-risk under always treat from gform_2 / g-risk under always treat from gform_1
RR1<- gform_2[["result"]][["g-form risk"]][[94]]/gform_1[["result"]][["g-form risk"]][[94]] #94 (47+47) refers to g-risk at time=46 
RR1 <-   0.5251526 /0.7216449   
round(RR1,2) # 0.73


# Estimate cumulative Hazard Ratio (cHR) as the ratio of log survivals in the two treatment strategies
# cHR= log(1- g-risk under always treat from gform_2) / log(1- g-risk under always treat from gform_1)
# Refer to ESM 2 for more details about the corresponding equation
cHR1<- log(1-gform_2[["result"]][["g-form risk"]][[94]])/log(1-gform_1[["result"]][["g-form risk"]][[94]])
cHR1<- log(1-0.5251526)/log(1-0.7216449)
round(cHR1,2) # 0.58


# Estimate bootstrap 95 % CI
# For this analysis we set t0=46 not 48.
# It was not possible at time=48 months as no observations left afterwards 
bs_1<-gform_1$bootests
#dim(bs_1)
#head(bs_1)
#tail(bs_1)
#plot(bs_1[,2])
est_1<-bs_1[bs_1$t0==46,c(2,3)] 
est_1
summary(est_1$`Never treat`)
#head(est_1)
#dim(est_1)


bs_2<-gform_2$bootests
est_2<-bs_2[bs_2$t0==46,c(2,3)]
est_2
summary(est_2$`Always treat`) 


rr_est<-as.numeric(est_2[[1]]/est_1[[1]])
rr_est

chr_est<- as.numeric(log(1-est_2[[1]])/log(1-est_1[[1]]))
chr_est

quantile(rr_est,c(0.025,0.5,0.95), na.rm=T) # 0.6556173, 0.7885677
quantile(chr_est,c(0.025,0.5,0.95), na.rm=T)# 0.4846768,  0.6655341

mean(rr_est)  # 0.7263474
mean(chr_est) # 0.5826567
#
#
# 
# Specification 2: Apply gformula in Control and Exp arm with Restricted models----
#
#
#
# Specification 2: Apply gformula in Control with restricted models

df_gformula<- df_ipcw_pp
final.gformula_control<- subset(df_gformula, df_gformula$plannedtrt==0)
final.gformula_control<- as.data.table(final.gformula_control)
#
time_since_switch<- function(pool, histvars, time_name, t,id_name){
  if(t== 0){
    pool[progtd==1 & get(time_name)==t,':=' (progtime=0, now.progtd=1)]
    pool[progtd==0 & get(time_name)==t,':=' (progtime=NA, now.progtd=NA)]
  } else{
    pool[progtd==1 & lag1_progtd==0 & get(time_name)==t,':=' (progtime=t, now.progtd=1)]
    pool[progtd==1 & lag1_progtd==1 & get(time_name)==t,':=' (now.progtd=0)]
  }
}
#
obs_data = final.gformula_control
#
outcome_name = 'deathtd'
outcome_type = 'survival'
censor_name= 'censoringtd_new'
id='id'
time_points= 47
time_name= 'time'
covnames = c('progtd','ecogtdGR1','icprogtd','tltd_wins','trttd')
covtypes = c('absorbing','binary','absorbing','normal','absorbing')
histories = c(lagged,time_since_switch)
histvars = list(c('progtd','ecogtdGR1','icprogtd','tltd_wins','trttd'),c('progtd'))
#
covparams= list(covmodels= c(
  progtd ~ time  + RACEGR1 + ECOGGR1  + lag1_icprogtd + lag1_tltd_wins, 
  
  ecogtdGR1 ~ time + AGE + I(AGE^2) + RACEGR1 + ECOGGR1 + SMKHISGR1+
    MICNSFL + PRRADYN + STRATRAN  + lag1_tltd_wins  + lag1_icprogtd + lag1_progtd + lag1_trttd,
  
  icprogtd ~ time + AGE + RACEGR1 + PRRADYN + STRATRAN + lag1_progtd + progtd,
  
  tltd_wins ~  time + I(time^2) + AGE  + SEX + RACEGR1 + ECOGGR1 + MICNSFL + PRRADYN + STRATRAN +
    progtd  + lag1_trttd,
  
  trttd ~  time + AGE + I(AGE^2) + RACEGR1 + STRATRAN  + lag1_icprogtd +  progtime + I(progtime^2))) 

ymodel <- deathtd ~ time + I(time^2) + AGE + I(AGE^2) + SMKHISGR1  + STRATRAN +
  lag1_trttd + trttd + progtime

censor_model<- censoringtd_new ~ time + I(time^2) + AGE +  I(AGE^2)+ SEX + RACEGR1+
  ECOGGR1 + SMKHISGR1  + PRRADYN + STRATRAN + progtd+ lag1_trttd + trttd

intvars= list('trttd')
interventions <- list(list(c(static, rep(0,time_points))))
int_descript <- c('Never treat')
nsimul <- 1100 
nsamples<- 20 #500
ncores <- parallel::detectCores()-1
#
gform_3 <- gformula(obs_data = obs_data,
                    outcome_type = outcome_type,
                    id=id,
                    time_points = time_points,
                    time_name = time_name,
                    covnames = covnames,
                    outcome_name = outcome_name,
                    covtypes = covtypes, covparams = covparams,
                    ymodel = ymodel, intvars = intvars,
                    interventions = interventions,
                    int_descript = int_descript,histories = histories,
                    histvars = histvars,
                    nsimul = nsimul,
                    seed = 1234,
                    nsamples= nsamples,# 20
                    parallel=TRUE,
                    ncores=ncores,
                    ci_method = 'percentile',
                    boot_diag = TRUE,
                    model_fits = TRUE,
                    censor_name=censor_name,
                    censor_model=censor_model,
                    ipw_cutoff_quantile = 0.95,
                    basecovs = c('AGE','SEX','RACEGR1','ECOGGR1','SMKHISGR1','MICNSFL','PRRADYN','STRATRAN'),
                    restrictions = list(c('trttd','progtd==1',simple_restriction,0)),
                    yrestrictions = list(c('progtd==1',0)),
                    sim_data_b = F)
print(gform_3)
#
#
#
# Specification 2: Apply gformula in Exp arm with restricted models

df_gformula<- df_ipcw_pp
final.gformula_exp<- subset(df_gformula, df_gformula$plannedtrt==1)
final.gformula_exp<- as.data.table(final.gformula_exp)
#
time_since_switch<- function(pool, histvars, time_name, t,id_name){
  if(t== 0){
    pool[progtd==1 & get(time_name)==t,':=' (progtime=0, now.progtd=1)]
    pool[progtd==0 & get(time_name)==t,':=' (progtime=NA, now.progtd=NA)]
  } else{
    pool[progtd==1 & lag1_progtd==0 & get(time_name)==t,':=' (progtime=t, now.progtd=1)]
    pool[progtd==1 & lag1_progtd==1 & get(time_name)==t,':=' (now.progtd=0)]
  }
}
#
obs_data = final.gformula_exp
class(obs_data)
outcome_name = 'deathtd'
outcome_type = 'survival'
censor_name= 'censoringtd_new'
id='id'
time_points= 47
time_name= 'time'
covnames = c('progtd','ecogtdGR1','icprogtd','tltd_wins','trttd')
covtypes = c('absorbing','binary','absorbing','normal','binary')
histories = c(lagged,time_since_switch)
histvars = list(c('progtd','ecogtdGR1','icprogtd','tltd_wins','trttd'),c('progtd'))
covparams= list(covmodels= c(progtd ~ time  + RACEGR1 + ECOGGR1 + SMKHISGR1 + STRATRAN +
                               lag1_icprogtd + lag1_tltd_wins,
                             
                             ecogtdGR1 ~ time + I(time^2) + AGE + I(AGE^2) +SEX + RACEGR1 + ECOGGR1 +SMKHISGR1+
                               MICNSFL + PRRADYN + STRATRAN  +
                               lag1_tltd_wins + lag1_icprogtd , 
                             
                             icprogtd ~ time  + ECOGGR1 + STRATRAN  + progtd +
                               lag1_progtd ,
                             
                             tltd_wins ~ time + I(time^2) + AGE +  I(AGE^2) + SEX 
                             + RACEGR1 + ECOGGR1 + STRATRAN + progtd ,
                             
                             trttd ~ 1)) 

ymodel <- deathtd ~  time  + SEX + MICNSFL + STRATRAN  +
  icprogtd   + lag1_icprogtd + progtime

censor_model<- censoringtd_new ~ time + I(time^2)+ ECOGGR1 + SMKHISGR1 + progtd + lag1_ecogtdGR1

intvars= list('trttd')
interventions <- list(list(c(static, rep(1,time_points))))
int_descript <- c('Always treat')
nsimul <- 1100 
nsamples<- 20 # 500
ncores <- parallel::detectCores()-1
#
gform_4 <- gformula(obs_data = obs_data,
                    outcome_type = outcome_type,
                    id=id,
                    time_points = time_points,
                    time_name = time_name,
                    covnames = covnames,
                    outcome_name = outcome_name,
                    covtypes = covtypes, 
                    covparams = covparams,
                    ymodel = ymodel,
                    intvars = intvars,
                    interventions = interventions,
                    int_descript = int_descript,
                    histories = histories,
                    histvars = histvars,
                    basecovs = c('AGE','SEX','RACEGR1','ECOGGR1','SMKHISGR1','MICNSFL','PRRADYN','STRATRAN'),                    
                    yrestrictions = list(c('progtd==1',0)),
                    nsimul = nsimul,
                    seed = 1234,
                    nsamples= nsamples,
                    parallel=TRUE,
                    ncores=ncores,
                    ci_method = 'percentile',
                    censor_name=censor_name,
                    censor_model=censor_model,
                    ipw_cutoff_quantile = 0.95,#
                    boot_diag = TRUE,
                    model_fits = TRUE)
print(gform_4)
#
#
#
#### Specification 2: RR and cHR------------------------------------------------

RR2<- gform_4[["result"]][["g-form risk"]][[94]]/gform_3[["result"]][["g-form risk"]][[94]] # 94 (47+47) refers to g-risk at time=46 
RR2 <-  0.5388889 /0.7114387        
round(RR2,2)  # 0.76


cHR2<- log(1-gform_4[["result"]][["g-form risk"]][[94]])/log(1-gform_3[["result"]][["g-form risk"]][[94]])
cHR2 <-  log(1-0.5388889)/log(1-0.7114387 )
round(cHR2,2) #  0.62


#Bootstrap CI
#
bs_1<-gform_3$bootests
est_1<-bs_1[bs_1$t0==46,c(2,3)]
est_1
summary(est_1$`Never treat`)


bs_2<-gform_4$bootests
est_2<-bs_2[bs_2$t0==46,c(2,3)]
est_2
summary(est_2$`Always treat`)


rr_est<-as.numeric(est_2[[1]]/est_1[[1]])
rr_est

chr_est<- as.numeric(log(1-est_2[[1]])/log(1-est_1[[1]]))
chr_est
#
quantile(rr_est,c(0.025,0.5,0.95),na.rm = T) # 0.6780236, 0.8143614 
quantile(chr_est,c(0.025,0.5,0.95),na.rm = T)# 0.5263397, 0.6988761 
#
mean(rr_est) #  0.7451444
mean(chr_est)#  0.613767

# Specification 3: Apply gformula in Control and Exp arm with further Restricted models-------
#
#
#
# Specification 3: Apply gformula in Control with further Restricted models

df_gformula<- df_ipcw_pp
final.gformula_control<- subset(df_gformula, df_gformula$plannedtrt==0)
final.gformula_control<- as.data.table(final.gformula_control)
#
time_since_switch<- function(pool, histvars, time_name, t,id_name){
  if(t== 0){
    pool[progtd==1 & get(time_name)==t,':=' (progtime=0, now.progtd=1)]
    pool[progtd==0 & get(time_name)==t,':=' (progtime=NA, now.progtd=NA)]
  } else{
    pool[progtd==1 & lag1_progtd==0 & get(time_name)==t,':=' (progtime=t, now.progtd=1)]
    pool[progtd==1 & lag1_progtd==1 & get(time_name)==t,':=' (now.progtd=0)]
  }
}
#
obs_data = final.gformula_control
#
outcome_name = 'deathtd'
outcome_type = 'survival'
censor_name= 'censoringtd_new'
id='id'
time_points= 47
time_name= 'time'
covnames = c('progtd','ecogtdGR1','icprogtd','tltd_wins','trttd')
covtypes = c('absorbing','binary','absorbing','normal','absorbing')
histories = c(lagged,time_since_switch)
histvars = list(c('progtd','ecogtdGR1','icprogtd','tltd_wins','trttd'),c('progtd'))
#
covparams= list(covmodels= c(
  progtd ~ time  + RACEGR1 + ECOGGR1  + lag1_icprogtd + lag1_tltd_wins, 
  
  ecogtdGR1 ~ time + AGE + I(AGE^2) + RACEGR1 + ECOGGR1 + SMKHISGR1+
    MICNSFL + PRRADYN + STRATRAN  + lag1_tltd_wins  + lag1_icprogtd + lag1_progtd + lag1_trttd,
  
  icprogtd ~ time + AGE  + PRRADYN + STRATRAN + lag1_progtd + progtd,
  
  tltd_wins ~  time + I(time^2) + AGE  + SEX + RACEGR1 + ECOGGR1 + MICNSFL + PRRADYN + STRATRAN +
    progtd  + lag1_trttd,
  
  trttd ~  time + RACEGR1 + STRATRAN  + lag1_icprogtd +  progtime)) 

ymodel <- deathtd ~ time  + AGE + I(AGE^2)  + STRATRAN +
  lag1_trttd + trttd + progtime

censor_model<- censoringtd_new ~ time + I(time^2) + AGE +  I(AGE^2)+ SEX + RACEGR1+ECOGGR1 +
  SMKHISGR1  + PRRADYN + progtd + lag1_trttd + trttd

intvars= list('trttd')
interventions <- list(list(c(static, rep(0,time_points))))
int_descript <- c('Never treat')
nsimul <- 1100 
nsamples<- 20 # 500
ncores <- parallel::detectCores()-1
#
gform_5 <- gformula(obs_data = obs_data,
                    outcome_type = outcome_type,
                    id=id,
                    time_points = time_points,
                    time_name = time_name,
                    covnames = covnames,
                    outcome_name = outcome_name,
                    covtypes = covtypes, covparams = covparams,
                    ymodel = ymodel, intvars = intvars,
                    interventions = interventions,
                    int_descript = int_descript,histories = histories,
                    histvars = histvars,
                    nsimul = nsimul,
                    seed = 1234,
                    nsamples= nsamples,# 20
                    parallel=TRUE,
                    ncores=ncores,
                    ci_method = 'percentile',
                    boot_diag = TRUE,
                    model_fits = TRUE,
                    censor_name=censor_name,
                    censor_model=censor_model,
                    ipw_cutoff_quantile = 0.95,
                    basecovs = c('AGE','SEX','RACEGR1','ECOGGR1','SMKHISGR1','MICNSFL','PRRADYN','STRATRAN'),
                    restrictions = list(c('trttd','progtd==1',simple_restriction,0)),
                    yrestrictions = list(c('progtd==1',0)),
                    sim_data_b = F)
print(gform_5)
#
#
#
# Specification 3: Apply gformula in Exp arm with further restricted models

df_gformula<- df_ipcw_pp
final.gformula_exp<- subset(df_gformula, df_gformula$plannedtrt==1)
final.gformula_exp<- as.data.table(final.gformula_exp)
#
time_since_switch<- function(pool, histvars, time_name, t,id_name){
  if(t== 0){
    pool[progtd==1 & get(time_name)==t,':=' (progtime=0, now.progtd=1)]
    pool[progtd==0 & get(time_name)==t,':=' (progtime=NA, now.progtd=NA)]
  } else{
    pool[progtd==1 & lag1_progtd==0 & get(time_name)==t,':=' (progtime=t, now.progtd=1)]
    pool[progtd==1 & lag1_progtd==1 & get(time_name)==t,':=' (now.progtd=0)]
  }
}
#
obs_data = final.gformula_exp
class(obs_data)
outcome_name = 'deathtd'
outcome_type = 'survival'
censor_name= 'censoringtd_new'
id='id'
time_points= 47
time_name= 'time'
covnames = c('progtd','ecogtdGR1','icprogtd','tltd_wins','trttd')
covtypes = c('absorbing','binary','absorbing','normal','binary')
histories = c(lagged,time_since_switch)
histvars = list(c('progtd','ecogtdGR1','icprogtd','tltd_wins','trttd'),c('progtd'))
covparams= list(covmodels= c(progtd ~ time  + RACEGR1 + ECOGGR1 + STRATRAN +
                               lag1_icprogtd + lag1_tltd_wins,
                             
                            ecogtdGR1 ~ time + AGE + I(AGE^2) +SEX + RACEGR1 + ECOGGR1 +SMKHISGR1+
                               MICNSFL + PRRADYN + STRATRAN + lag1_icprogtd , 
                             
                             icprogtd ~ time  + ECOGGR1 + STRATRAN  + progtd +
                               lag1_progtd ,
                             
                             tltd_wins ~ time + I(time^2) + AGE +  I(AGE^2) + SEX 
                             + RACEGR1 + ECOGGR1 + STRATRAN +  + progtd ,
                             
                             trttd ~ 1)) 
ymodel <- deathtd ~  time   + icprogtd   + lag1_icprogtd + progtime

censor_model<- censoringtd_new ~ time + I(time^2)+ ECOGGR1 + SMKHISGR1 + progtd + lag1_ecogtdGR1

intvars= list('trttd')
interventions <- list(list(c(static, rep(1,time_points))))
int_descript <- c('Always treat')
nsimul <- 1100
nsamples<- 20 # 500
ncores <- parallel::detectCores()-1
#
gform_6 <- gformula(obs_data = obs_data,
                    outcome_type = outcome_type,
                    id=id,
                    time_points = time_points,
                    time_name = time_name,
                    covnames = covnames,
                    outcome_name = outcome_name,
                    covtypes = covtypes, 
                    covparams = covparams,
                    ymodel = ymodel,
                    intvars = intvars,
                    interventions = interventions,
                    int_descript = int_descript,
                    histories = histories,
                    histvars = histvars,
                    basecovs = c('AGE','SEX','RACEGR1','ECOGGR1','SMKHISGR1','MICNSFL','PRRADYN','STRATRAN'),               
                    yrestrictions = list(c('progtd==1',0)),
                    nsimul = nsimul,
                    seed = 1234,
                    nsamples= nsamples,
                    parallel=TRUE,
                    ncores=ncores,
                    ci_method = 'percentile',
                    censor_name=censor_name,
                    censor_model=censor_model,
                    ipw_cutoff_quantile = 0.95,#
                    boot_diag = TRUE,
                    model_fits = TRUE)
print(gform_6)
# 
#
#
#### Specification 3: RR and cHR------------------------------------------------

RR3<- gform_6[["result"]][["g-form risk"]][[94]]/gform_5[["result"]][["g-form risk"]][[94]] # 94 (47+47) refers to g-risk at time=46 
RR3 <-  0.5240447/0.7111855           
round(RR3,2) # 0.74


cHR3<- log(1-gform_6[["result"]][["g-form risk"]][[94]])/log(1-gform_5[["result"]][["g-form risk"]][[94]])
cHR3 <- log(1-0.5240447)/log(1-0.7111855)
round(cHR3,2) # 0.60


# Bootstrap CI
bs_1<-gform_5$bootests
est_1<-bs_1[bs_1$t0==46,c(2,3)]
est_1
summary(est_1$`Never treat`)


bs_2<-gform_6$bootests
est_2<-bs_2[bs_2$t0==46,c(2,3)]
est_2
summary(est_2$`Always treat`)

rr_est<-as.numeric(est_2[[1]]/est_1[[1]])
rr_est

chr_est<- as.numeric(log(1-est_2[[1]])/log(1-est_1[[1]]))
chr_est

quantile(rr_est,c(0.025,0.5,0.95), na.rm = T) # 0.6775300, 0.8018994 
quantile(chr_est,c(0.025,0.5,0.95), na.rm = T)# 0.5274721, 0.6830030
#
mean(rr_est) # 0.7469178
mean(chr_est)# 0.6143779
#
#
#
# Specification 4: sensitivity with ecogtd012-----------------------------------
# Same as spec 1 but replace binary ecogtGR1 with categorical ecogtdGR2
#
#
#
# Specification 4: sensitivity with ecogtd012 in control

df_gformula<- df_ipcw_pp
final.gformula_control<- subset(df_gformula, df_gformula$plannedtrt==0)
final.gformula_control<- as.data.table(final.gformula_control)
#
time_since_switch<- function(pool, histvars, time_name, t,id_name){
  if(t== 0){
    pool[progtd==1 & get(time_name)==t,':=' (progtime=0, now.progtd=1)]
    pool[progtd==0 & get(time_name)==t,':=' (progtime=NA, now.progtd=NA)]
  } else{
    pool[progtd==1 & lag1_progtd==0 & get(time_name)==t,':=' (progtime=t, now.progtd=1)]
    pool[progtd==1 & lag1_progtd==1 & get(time_name)==t,':=' (now.progtd=0)]
  }
}
#
obs_data = final.gformula_control
#
outcome_name = 'deathtd'
outcome_type = 'survival'
censor_name= 'censoringtd_new'
id='id'
time_points= 47
time_name= 'time'
covnames = c('progtd','ecogtdGR2','icprogtd','tltd_wins','trttd')
covtypes = c('absorbing','categorical','absorbing','normal','absorbing')
histories = c(lagged,time_since_switch)
histvars = list(c('progtd','ecogtdGR2','icprogtd','tltd_wins','trttd'),c('progtd'))

covparams= list(covmodels= c(
  progtd ~ time + I(time^2) + RACEGR1 + ECOGGR1 + SMKHISGR1 + lag1_icprogtd + lag1_tltd_wins, 
  
  ecogtdGR2 ~ time + AGE + I(AGE^2) + SEX + RACEGR1 + ECOGGR1 + SMKHISGR1+
    MICNSFL + PRRADYN + STRATRAN  + lag1_tltd_wins  + lag1_icprogtd + lag1_progtd + progtd + lag1_trttd,
  
  icprogtd ~ time + AGE + RACEGR1 + PRRADYN + STRATRAN + lag1_progtd + progtd + lag1_tltd_wins + lag1_trttd,
  
  tltd_wins ~  time + I(time^2) + AGE  + SEX + RACEGR1 + ECOGGR1 + MICNSFL + PRRADYN + STRATRAN +
    lag1_ecogtdGR2 +lag1_icprogtd +lag1_progtd +ecogtdGR2 + progtd + icprogtd + lag1_trttd,
  
  trttd ~  time + AGE + I(AGE^2) + RACEGR1 + STRATRAN  + lag1_icprogtd + icprogtd +  progtime + I(progtime^2))) 

ymodel <- deathtd ~ time + I(time^2) + AGE + I(AGE^2) + SMKHISGR1  + STRATRAN  + lag1_ecogtdGR2 +
  ecogtdGR2  + lag1_trttd + trttd + progtime

censor_model<- censoringtd_new ~ time + I(time^2) + AGE +  I(AGE^2)+ SEX + RACEGR1+
  ECOGGR1 + SMKHISGR1  + PRRADYN + STRATRAN  + lag1_icprogtd + icprogtd +
  lag1_tltd_wins + tltd_wins + lag1_progtd + progtd + ecogtdGR2 + lag1_ecogtdGR2 +
  lag1_trttd + trttd # + progtime (not possible, because of NA)

intvars= list('trttd')
interventions <- list(list(c(static, rep(0,time_points))))
int_descript <- c('Never treat')
nsimul <- 1100 
nsamples<- 20   # 500
ncores <- parallel::detectCores()-1
#
gform_7 <- gformula(obs_data = obs_data,
                    outcome_type = outcome_type,
                    id=id,
                    time_points = time_points,
                    time_name = time_name,
                    covnames = covnames,
                    outcome_name = outcome_name,
                    covtypes = covtypes, covparams = covparams,
                    ymodel = ymodel, intvars = intvars,
                    interventions = interventions,
                    int_descript = int_descript,histories = histories,
                    histvars = histvars,
                    nsimul = nsimul,
                    seed = 1234,
                    nsamples= nsamples,# 20
                    parallel=TRUE,
                    ncores=ncores,
                    ci_method = 'percentile',
                    boot_diag = TRUE,
                    model_fits = TRUE,
                    censor_name=censor_name,
                    censor_model=censor_model,
                    ipw_cutoff_quantile = 0.95,
                    basecovs = c('AGE','SEX','RACEGR1','ECOGGR1','SMKHISGR1','MICNSFL','PRRADYN','STRATRAN'),
                    restrictions = list(c('trttd','progtd==1',simple_restriction,0)),
                    yrestrictions = list(c('progtd==1',0)),
                    sim_data_b = F)
print(gform_7)
#
#
#
# Specification 4: sensitivity with ecogtd012 in exp

df_gformula<- df_ipcw_pp
final.gformula_exp<- subset(df_gformula, df_gformula$plannedtrt==1)
final.gformula_exp<- as.data.table(final.gformula_exp)
#
time_since_switch<- function(pool, histvars, time_name, t,id_name){
  if(t== 0){
    pool[progtd==1 & get(time_name)==t,':=' (progtime=0, now.progtd=1)]
    pool[progtd==0 & get(time_name)==t,':=' (progtime=NA, now.progtd=NA)]
  } else{
    pool[progtd==1 & lag1_progtd==0 & get(time_name)==t,':=' (progtime=t, now.progtd=1)]
    pool[progtd==1 & lag1_progtd==1 & get(time_name)==t,':=' (now.progtd=0)]
  }
}
#
obs_data = final.gformula_exp
class(obs_data)
outcome_name = 'deathtd'
outcome_type = 'survival'
censor_name= 'censoringtd_new'
id='id'
time_points= 47
time_name= 'time'
covnames = c('progtd','ecogtdGR2','icprogtd','tltd_wins','trttd')
covtypes = c('absorbing','categorical','absorbing','normal','binary')
histories = c(lagged,time_since_switch)
histvars = list(c('progtd','ecogtdGR2','icprogtd','tltd_wins','trttd'),c('progtd'))

covparams= list(covmodels= c(progtd ~ time + SEX + RACEGR1 + ECOGGR1 + SMKHISGR1 + STRATRAN +
                               lag1_icprogtd + lag1_tltd_wins,
                             
                             ecogtdGR2 ~ time + I(time^2) + AGE + I(AGE^2) +SEX + RACEGR1 + ECOGGR1 +SMKHISGR1+
                               MICNSFL + PRRADYN + STRATRAN  +
                               lag1_tltd_wins + lag1_icprogtd , 
                             
                             icprogtd ~ time + I(time^2) + ECOGGR1 + STRATRAN  + progtd +
                               lag1_progtd + lag1_tltd_wins ,
                             
                             tltd_wins ~ time + I(time^2) + AGE +  I(AGE^2) + SEX 
                             + RACEGR1 + ECOGGR1 + STRATRAN + ecogtdGR2 + progtd + lag1_ecogtdGR2 +
                               lag1_progtd + lag1_icprogtd + icprogtd,
                             
                             trttd ~ 1)) 

ymodel <- deathtd ~  time  + SEX + MICNSFL + STRATRAN + ecogtdGR2 +
  icprogtd  +lag1_ecogtdGR2 + lag1_icprogtd + progtime

censor_model<- censoringtd_new ~ time + I(time^2)+ ECOGGR1 + SMKHISGR1 +
  ecogtdGR2 + icprogtd  + progtd + lag1_ecogtdGR2 + lag1_icprogtd + lag1_progtd

intvars= list('trttd')
interventions <- list(list(c(static, rep(1,time_points))))
int_descript <- c('Always treat')
nsimul <- 1100
nsamples<- 20 # 500 
ncores <- parallel::detectCores()-1
#
gform_8 <- gformula(obs_data = obs_data,
                    outcome_type = outcome_type,
                    id=id,
                    time_points = time_points,
                    time_name = time_name,
                    covnames = covnames,
                    outcome_name = outcome_name,
                    covtypes = covtypes, 
                    covparams = covparams,
                    ymodel = ymodel,
                    intvars = intvars,
                    interventions = interventions,
                    int_descript = int_descript,
                    histories = histories,
                    histvars = histvars,
                    basecovs = c('AGE','SEX','RACEGR1','ECOGGR1','SMKHISGR1','MICNSFL','PRRADYN','STRATRAN'),                   
                    yrestrictions = list(c('progtd==1',0)),
                    nsimul = nsimul,
                    seed = 1234,
                    nsamples= nsamples,
                    parallel=TRUE,
                    ncores=ncores,
                    ci_method = 'percentile',
                    censor_name=censor_name,
                    censor_model=censor_model,
                    ipw_cutoff_quantile = 0.95,#
                    boot_diag = TRUE,
                    model_fits = TRUE)
print(gform_8) 
#
#
#
#### Specification 4: RR and cHR------------------------------------------------

RR4<- gform_8[["result"]][["g-form risk"]][[94]]/gform_7[["result"]][["g-form risk"]][[94]] # 94 (47+47) refers to g-risk at time=46 
RR4 <-  0.5200147  /0.7261180        
round(RR4,2) # 0.72


cHR4<- log(1-gform_8[["result"]][["g-form risk"]][[94]])/log(1-gform_7[["result"]][["g-form risk"]][[94]])
cHR4 <-  log(1-0.5200147 )/log(1-0.7261180) 
round(cHR4,2) # 0.57


#Bootstrap CI
#
bs_1<-gform_7$bootests
est_1<-bs_1[bs_1$t0==46,c(2,3)]
est_1
summary(est_1$`Never treat`)


bs_2<-gform_8$bootests
est_2<-bs_2[bs_2$t0==46,c(2,3)]
est_2
summary(est_2$`Always treat`)# 61 NA


rr_est<-as.numeric(est_2[[1]]/est_1[[1]])
rr_est

chr_est<- as.numeric(log(1-est_2[[1]])/log(1-est_1[[1]]))
chr_est

quantile(rr_est,c(0.025,0.5,0.95),na.rm = T)# 0.6408039, 0.7742088 
quantile(chr_est,c(0.025,0.5,0.95),na.rm = T)# 0.4770717,0.6539063
#
mean(rr_est)# 0.7168801
mean(chr_est)# 0.5765526
#
#
#
# Specification 5: sensitivity without progtime---------------------------------
# Same as spec 1 but remove progtime from the outcome and trt model
#
#
#
# Specification 5: sensitivity without progtime in control

df_gformula<- df_ipcw_pp
final.gformula_control<- subset(df_gformula, df_gformula$plannedtrt==0)
final.gformula_control<- as.data.table(final.gformula_control)
#
obs_data = final.gformula_control
summary(obs_data)
class(obs_data)
#
outcome_name = 'deathtd'
outcome_type = 'survival'
censor_name= 'censoringtd_new'
id='id'
time_points= 47
time_name= 'time'
covnames = c('progtd','ecogtdGR1','icprogtd','tltd_wins','trttd')
covtypes = c('absorbing','binary','absorbing','normal','absorbing')
histories = c(lagged)
histvars = list(c('progtd','ecogtdGR1','icprogtd','tltd_wins','trttd'))
covparams= list(covmodels= c(
  progtd ~ time + I(time^2) + RACEGR1 + ECOGGR1 + SMKHISGR1 + lag1_icprogtd + lag1_tltd_wins, 
  
  ecogtdGR1 ~ time + AGE + I(AGE^2) + SEX + RACEGR1 + ECOGGR1 + SMKHISGR1+
    MICNSFL + PRRADYN + STRATRAN  + lag1_tltd_wins  + lag1_icprogtd + lag1_progtd + progtd + lag1_trttd,
  
  icprogtd ~ time + AGE + RACEGR1 + PRRADYN + STRATRAN + lag1_progtd + progtd + lag1_tltd_wins + lag1_trttd,
  
  tltd_wins ~  time + I(time^2) + AGE  + SEX + RACEGR1 + ECOGGR1 + MICNSFL + PRRADYN + STRATRAN +
    lag1_ecogtdGR1 +lag1_icprogtd +lag1_progtd +ecogtdGR1 + progtd + icprogtd + lag1_trttd,
  
  trttd ~  time + AGE + I(AGE^2) + RACEGR1 + STRATRAN  + lag1_icprogtd + icprogtd)) 

ymodel <- deathtd ~ time + I(time^2) + AGE + I(AGE^2) + SMKHISGR1  + STRATRAN  + lag1_ecogtdGR1 +
  ecogtdGR1  + lag1_trttd + trttd 

censor_model<- censoringtd_new ~ time + I(time^2) + AGE +  I(AGE^2)+ SEX + RACEGR1+
  ECOGGR1 + SMKHISGR1  + PRRADYN + STRATRAN  + lag1_icprogtd + icprogtd +
  lag1_tltd_wins + tltd_wins + lag1_progtd + progtd + ecogtdGR1 + lag1_ecogtdGR1 +
  lag1_trttd + trttd 

intvars= list('trttd')
interventions <- list(list(c(static, rep(0,time_points))))
int_descript <- c('Never treat')
nsimul <- 1100
nsamples<- 20 # 500   
ncores <- parallel::detectCores()-1
#
gform_9<- gformula(obs_data = obs_data,
                   outcome_type = outcome_type,
                   id=id,
                   time_points = time_points,
                   time_name = time_name,
                   covnames = covnames,
                   outcome_name = outcome_name,
                   covtypes = covtypes, covparams = covparams,
                   ymodel = ymodel, intvars = intvars,
                   interventions = interventions,
                   int_descript = int_descript,histories = histories,
                   histvars = histvars,
                   nsimul = nsimul,
                   seed = 1234,
                   nsamples= nsamples,# 20
                   parallel=TRUE,
                   ncores=ncores,
                   ci_method = 'percentile',
                   boot_diag = TRUE,
                   model_fits = TRUE,
                   censor_name=censor_name,
                   censor_model=censor_model,
                   ipw_cutoff_quantile = 0.95,
                   basecovs = c('AGE','SEX','RACEGR1','ECOGGR1','SMKHISGR1','MICNSFL','PRRADYN','STRATRAN'),
                   restrictions = list(c('trttd','progtd==1',simple_restriction,0)),
                   yrestrictions = list(c('progtd==1',0)),
                   sim_data_b = F)
print(gform_9) 
# 
#
#
# Specification 5: sensitivity without progtime in exp

df_gformula<- df_ipcw_pp
final.gformula_exp<- subset(df_gformula, df_gformula$plannedtrt==1)
final.gformula_exp<- as.data.table(final.gformula_exp)
#
time_since_switch<- function(pool, histvars, time_name, t,id_name){
  if(t== 0){
    pool[progtd==1 & get(time_name)==t,':=' (progtime=0, now.progtd=1)]
    pool[progtd==0 & get(time_name)==t,':=' (progtime=NA, now.progtd=NA)]
  } else{
    pool[progtd==1 & lag1_progtd==0 & get(time_name)==t,':=' (progtime=t, now.progtd=1)]
    pool[progtd==1 & lag1_progtd==1 & get(time_name)==t,':=' (now.progtd=0)]
  }
}
#
obs_data = final.gformula_exp
class(obs_data)
outcome_name = 'deathtd'
outcome_type = 'survival'
censor_name= 'censoringtd_new'
id='id'
time_points= 47
time_name= 'time'
covnames = c('progtd','ecogtdGR1','icprogtd','tltd_wins','trttd')
covtypes = c('absorbing','binary','absorbing','normal','binary')
histories = c(lagged,time_since_switch)
histvars = list(c('progtd','ecogtdGR1','icprogtd','tltd_wins','trttd'),c('progtd'))
covparams= list(covmodels= c(progtd ~ time + SEX + RACEGR1 + ECOGGR1 + SMKHISGR1 + STRATRAN +
                               lag1_icprogtd + lag1_tltd_wins,
                             
                             ecogtdGR1 ~ time + I(time^2) + AGE + I(AGE^2) +SEX + RACEGR1 + ECOGGR1 +SMKHISGR1+
                               MICNSFL + PRRADYN + STRATRAN  +
                               lag1_tltd_wins + lag1_icprogtd , 
                             
                             icprogtd ~ time + I(time^2) + ECOGGR1 + STRATRAN  + progtd +
                               lag1_progtd + lag1_tltd_wins ,
                             
                             tltd_wins ~ time + I(time^2) + AGE +  I(AGE^2) + SEX 
                             + RACEGR1 + ECOGGR1 + STRATRAN + ecogtdGR1 + progtd + lag1_ecogtdGR1 +
                               lag1_progtd + lag1_icprogtd + icprogtd,
                             
                             trttd ~ 1)) 

ymodel <- deathtd ~  time  + SEX + MICNSFL + STRATRAN + ecogtdGR1 +
  icprogtd  +lag1_ecogtdGR1 + lag1_icprogtd 

censor_model<- censoringtd_new ~ time + I(time^2)+ ECOGGR1 + SMKHISGR1 +
  ecogtdGR1 + icprogtd  + progtd + lag1_ecogtdGR1 + lag1_icprogtd + lag1_progtd

intvars= list('trttd')
interventions <- list(list(c(static, rep(1,time_points))))
int_descript <- c('Always treat')
nsimul <- 1100 
nsamples<- 20 # 500
ncores <- parallel::detectCores()-1
#
gform_10<- gformula(obs_data = obs_data,
                   outcome_type = outcome_type,
                    id=id,
                    time_points = time_points,
                    time_name = time_name,
                    covnames = covnames,
                    outcome_name = outcome_name,
                    covtypes = covtypes, 
                    covparams = covparams,
                    ymodel = ymodel,
                    intvars = intvars,
                    interventions = interventions,
                    int_descript = int_descript,
                    histories = histories,
                    histvars = histvars,
                   basecovs = c('AGE','SEX','RACEGR1','ECOGGR1','SMKHISGR1','MICNSFL','PRRADYN','STRATRAN'),                 
                    yrestrictions = list(c('progtd==1',0)),
                    nsimul = nsimul,
                    seed = 1234,
                    nsamples= nsamples,
                    parallel=TRUE,
                    ncores=ncores,
                    ci_method = 'percentile',
                    censor_name=censor_name,
                    censor_model=censor_model,
                    ipw_cutoff_quantile = 0.95,#
                    boot_diag = TRUE,
                    model_fits = TRUE)

print(gform_10)
#
#
#
### Specification 5: RR and cHR-------------------------------------------------

RR5<- gform_10[["result"]][["g-form risk"]][[94]]/gform_9[["result"]][["g-form risk"]][[94]] # 94 (47+47) refers to g-risk at time=46 
RR5 <-  0.5276746 /0.7154957      
round(RR5,2) # 0.74


cHR5<- log(1-gform_10[["result"]][["g-form risk"]][[94]])/log(1-gform_9[["result"]][["g-form risk"]][[94]])
cHR5 <-  log(1-0.5276746) /log(1-0.7154957)
round(cHR5,2)# 0.60


#Bootstrap CI
#
bs_1<-gform_9$bootests
est_1<-bs_1[bs_1$t0==46,c(2,3)]
est_1
summary(est_1$`Never treat`)


bs_2<-gform_10$bootests
est_2<-bs_2[bs_2$t0==46,c(2,3)]
est_2
summary(est_2$`Always treat`)# 1NA

rr_est<-as.numeric(est_2[[1]]/est_1[[1]])
rr_est

chr_est<- as.numeric(log(1-est_2[[1]])/log(1-est_1[[1]]))
chr_est
#
quantile(rr_est,c(0.025,0.5,0.95),na.rm=T)# 0.6568967, 0.7726305 
quantile(chr_est,c(0.025,0.5,0.95),na.rm=T)# 0.4923546,0.6345899 
#
mean(rr_est)# 0.7258088
mean(chr_est)# 0.5795808
#
#
#
# Forest plot for g-formula-----------------------------------------------------

# Create Data for the forest plot
# Create a data frame with specification, mean, lower bounds, upper bounds,
# cumulative hazard ratios (cHR), and 95% confidence intervals (CI).
data_forest<- data.frame(
  Specification= c("1","2","3","4","5"),
  mean= c(0.58,0.62,0.60,0.57,0.60),
  lower=c(0.48,0.53,0.53,0.48,0.49),
  upper=c(0.67,0.70,0.68,0.65,0.63),
  cHR= c("0.58","0.62","0.60","0.57","0.60"),
  CI= c("0.48,0.67","0.53,0.70","0.53,0.68","0.48,0.65","0.49,0.63")
)

# Generate Forest Plot:
# Use the 'forestplot' function to generate a forest plot.
forest_all<- data_forest |>
  forestplot(labeltext= c(Specification,cHR,CI),
             clip= c(0.40,1.2),
             vertices=TRUE,
             xlog=TRUE,
             xlab= "cHR (95% CI)",
             boxsize=0.1,
             graph.pos=4,
             text_size=2,
             spacing=0,
             zero=1,
             #xticks=c(0.1,0.2),
             txt_gp=fpTxtGp(cex=0.95),
             axes=gpar(cex=0.7),
             #at= c(0.2,0.4,0.6),
             #labels= c("0.2","0.4","0.6"),
             #title= "Forest plot: Cumulative Hazard ratio with 95% Confidence interval of Overall Survival adjusted for Lost to follow-up and treatment switching investigated using parametric g-formula"
  )|>
  fp_set_style(box="darkblue",
               line="black",
               summary=NULL)|>
  fp_add_header(Specification= c("","Specification"),
                cHR= c("","cHR"),
                CI= c("","95% CI"))|>
  fp_set_zebra_style("#EFEFEF")

# Display the Forest Plot:
forest_all
#
#
#
# Survival Plot ----------------------------------------------------------------

k=46

# Select which analysis (specification) you want to plot
# here we plot results from specification 1:  gform_1 and gform_2
result_control <- gform_1[["result"]]
result_0 <- subset(result_control, result_control$Interv.==1) 
result_0 <- result_0[,c("k","Interv.","g-form risk","Risk SE")]
result_0$Interv.<- 0
result_0$surv <- 1-result_0$`g-form risk`
result_0$se_s <- (result_0$surv)*(result_0$`Risk SE`)
result_0$ll<- (result_0$surv)-1.96*(result_0$se_s)
result_0$ul<- (result_0$surv)+1.96*(result_0$se_s)

result_exp <- gform_2[["result"]]
result_1 <- subset(result_exp, result_exp$Interv.==1) 
result_1 <- result_1[,c("k","Interv.","g-form risk","Risk SE")]
result_1$surv <- 1-result_1$`g-form risk`
result_1$se_s <- (result_1$surv)*(result_1$`Risk SE`)
result_1$ll<- (result_1$surv)-1.96*(result_1$se_s)
result_1$ul<- (result_1$surv)+1.96*(result_1$se_s)

results<- rbind(result_0,result_1)
results$Interv.<- factor(results$Interv.)
results$Interv.<- ifelse(results$Interv.==1,"brigatinib","crizotinib")


# Survival probability table:
results_prob<- subset(results, results$k==0|results$k==12|results$k==24|results$k==36|
                        results$k==42|results$k==46)
results_prob<- results_prob[,c("k","Interv.","surv","ll","ul")]
results_prob<- as.data.frame(results_prob)
results_prob[,sapply(results_prob, is.numeric)]<- 
  round(results_prob[,sapply(results_prob, is.numeric)],2)
results_prob[,c("surv","ll","ul")]<- results_prob[,c("surv","ll","ul")]*100
write.csv(results_prob, file="D:/../name.csv",row.names = T)


# Plot survival(%) over time 
results$surv<- results$surv*100
ggplot_gform<- ggplot(results,aes(x=k,y=surv))+
  geom_line(aes(colour=Interv.))+
  geom_point(aes(colour=Interv.))+
  xlab("Time after randomization (months)")+
  scale_x_continuous(limits=c(0,48),breaks=seq(0,48,6))+
  ylab(" Overall Survival (%)")+
  scale_y_continuous(limits=c(0,100),breaks=seq(0,100, by=25))+
  #ggtitle("Parametric G-formula: Overall Survival by each follow-up time")+
  labs(colour="")+
  theme_bw()+
  theme(legend.position = c(0.2,0.6),
        panel.border = element_blank(),
        legend.background = element_rect(fill="transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color="black"))+
  scale_color_manual(values=c("#F95700FF","#00539CFF"))
ggplot_gform

# add text (cHR and RR) from specification 1 to the ggplot
ggplot_gform + annotate("text", x=36,y=85, label="cHR 0.58 (95% CI 0.48,0.67) \n RR 0.73 (95% CI 0.66,0.79)", size=unit(4,"pt"))
#
#
#