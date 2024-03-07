# Set working directory:
setwd("D:\\LMU--2020\\PhD\\Publication\\Documents_for_revision")
#setwd("C:\\Users\\mansmann\\Documents\\IBE_Projekte\\Promotionen\\Promotion_Amani_Altawil\\Synthetic_Data")
#setwd("D:\\LMU--2020\\PhD\\Rscripts\\Synthetic data_UM")
setwd("C:\\Users\\altawil\\Desktop\\myproject_github")
getwd()
dir()

# The data is in a long format: for each patient 50 time points will be available.
N.pat<-2000     # Total number of patients in the two-armed trial
n.times<-50     # Number of time points for each patient
id<-rep(1:N.pat,rep(n.times,N.pat))     # Generate patient IDs in a long format
time<-rep((0):49,N.pat)     # Create a sequence of time points (0 to 49) for each patient
#
#
# Baseline Variables
#
#
# Generate data for treatment randomized at baseline
trt<-rbinom(N.pat,1,0.4982)     # Generate treatment assignment (binary: 0 or 1) for each patient
plannedtrt<-rep(trt,rep(n.times,N.pat))     # Repeat treatment assignment for each time point
#
#
#
#
#
# Generate data for binary variables:
#
# ICNSMTFL - brain metastasis Y/N (1/0)
# PCHEMFL  - previous chemotherapy Y/N (1/0)
# MICNSFL  - measurable intracranial CNS disease Y/N (1/0)
# ALKFDA   - test for cancer mutation Y/N (1/0)
# DIAGSTE  - diagnosis stage at study entry Y/N (1/0)
# PRACTDS  - prior anti-cancer therapy Y/N (1/0)
# PRRADYN  - prior radiation therapy Y/N (1/0)
# PRRBRYN  - Prior radiation therapy to the brain Y/N (1/0)
# RACEGR1  - race asian versus non-asian (1/0)
# SEX      - Female vs Males (1/0)
#
# Specify the probabilities for each variable
probabilities <- c(ICNSMTFL = 0.2945, PCHEMFL = 0.2655, MICNSFL = 0.1491, 
                   ALKFDA = 0.8545,DIAGSTE = 0.07273, PRACTDS = 0.2945,
                   PRRADYN = 0.2655, PRRBRYN = 0.1345,RACEGR1 = 0.3927, SEX = 0.5455)

# Create an empty list to store the variables
variable_list <- list()

# Generate random binary data (0 or 1) based on the specified probability
for (var_name in names(probabilities)) {
  xxx <- rbinom(N.pat, 1, probabilities[var_name])
  # Replicate the generated data for each patient and time point in the data set
  variable_list[[var_name]] <- rep(xxx, rep(n.times, N.pat))
}
df <- as.data.frame(variable_list)
#
#
#
#
#
# Generate data for categorical variables:
#
# Specify the frequencies of categorical variable
# Convert frequencies to probabilities
# Sample randomly for each patient
# Replicate sampled probabilities for each patient and time point in the data set
# Create different categories when needed
#
#
# ECOG - 0,1,2
p.ECOG<-c(107,154,14)     
p.ECOG<-p.ECOG/sum(p.ECOG)     
xxx<-sample(0:2,N.pat,replace=T,p.ECOG)     
ECOG<-rep(xxx,rep(n.times,N.pat))      
# ECOGGR1: 0 vs 1+
ECOGGR1<- ifelse(ECOG<1,0,1)     
# ECOGGR2: 01 vs 2
ECOGGR2<- ifelse(ECOG<2,0,1)     
#
#
# SMKHIS - Smoking history Never,former, current (0,1,2)
p.SMKHIS<-c(159,106,10)
p.SMKHIS<-p.SMKHIS/sum(p.SMKHIS)
xxx<-sample(0:2,N.pat,replace=T,p.SMKHIS)
SMKHIS<-rep(xxx,rep(n.times,N.pat))
# SMKHISGR1: never vs ever
SMKHISGR1<- ifelse(SMKHIS<1,0,1)  
#
#
# REGION1N - Region (0,1,2)
p.REGION1N<-c(25,143,107)
p.REGION1N<-p.REGION1N/sum(p.REGION1N)
xxx<-sample(0:2,N.pat,replace=T,p.REGION1N)
REGION1N<-rep(xxx,rep(n.times,N.pat))
#
#
# HISCLSE - tumor histopathologic class (Adenocarcinoma, Adenosquamouscarcinoma,largecell, other,squamous)
# coded as 0,1,2,3,4
p.HISCLSE<-c(263,4,2,2,4)
p.HISCLSE<-p.HISCLSE/sum(p.HISCLSE)
xxx<-sample(0:4,N.pat,replace=T,p.HISCLSE)
HISCLSE<-rep(xxx,rep(n.times,N.pat))
# HISCLSEGR1
HISCLSEGR1<- ifelse(HISCLSE==0|HISCLSE==2,0,1)     
# HISCLSEGR2
HISCLSEGR2<- factor(ifelse(HISCLSE==0|HISCLSE==2,0,
                    ifelse(HISCLSE==1|HISCLSE==4,1,2)))   
#
#
# DIAGINIT - initial diagnosis stage (IV,IIIB,IIIA,IIA,IB,IA) coded as 0,1,2,3,4,5
p.DIAGINIT<-c(213,21,21,7,7,6)
p.DIAGINIT<-p.DIAGINIT/sum(p.DIAGINIT)
xxx<-sample(0:5,N.pat,replace=T,p.DIAGINIT)
DIAGINIT<-rep(xxx,rep(n.times,N.pat))
# DIAGINITGR1
DIAGINITGR1<- factor(ifelse(DIAGINIT==0,0,
                            ifelse(DIAGINIT==1,1,2)))     
# DIAGINITGR2
DIAGINITGR2<- factor(ifelse(DIAGINIT==0,0,
                            ifelse(DIAGINIT==1|DIAGINIT==2,1,
                            ifelse(DIAGINIT==3,2,3))))   
# DIAGINITGR3
DIAGINITGR3<- factor(ifelse(DIAGINIT==0,0,
                            ifelse(DIAGINIT==1,1,
                            ifelse(DIAGINIT==2,2,3))))   
#
#
# STRATRAN: strata at randomization
# No iCNS Metastasis at baseline/No prior chemotherapy (0)
# iCNS Metastasis at baseline/No prior chemotherapy (1)
# No iCNS Metastasis at baseline/Prior chemotherapy (2)
# iCNS Metastasis at baseline/Prior chemotherapy (3)
p.STRATRAN<-c(149,53,45,28)
p.STRATRAN<-p.STRATRAN/sum(p.STRATRAN)
xxx<-sample(0:3,N.pat,replace=T,p.STRATRAN)
STRATRAN<-rep(xxx,rep(n.times,N.pat))
STRATRAN <- factor(STRATRAN)
#
#
# LIASE - Lung involvement at study entry (Both,Right,Left,Lung not involved) coded as 0,1,2,3
p.LIASE<-c(97,83,73,22)
p.LIASE<-p.LIASE/sum(p.LIASE)
xxx<-sample(0:3,N.pat,replace=T,p.LIASE)
LIASE<-rep(xxx,rep(n.times,N.pat))
LIASE <- factor(LIASE)
#
#
#
#
#
# Generate data for continuous variables:
#
# Specify the mean, SD, minimum and max of each variable
# Generate a random normal distribution of each variable
# Truncate the generated values to fall within the specified range [min, max]
# Repeat the values for each patient and time point in the data set
# Categorize some variables when needed
#
#
# AGE - Age
age.mean<-58.24 
age.sd<-12.03 
min.age <- 27     
max.age <- 89     
xxx<-rnorm(N.pat,age.mean,age.sd)    
xxx <- pmax(pmin(xxx, max.age), min.age)     
AGE<-rep(xxx,rep(n.times,N.pat))     
# Age 01
AGEGR1N<-ifelse(AGE>=65,1,0)     
#
#
# tl_b_wins - target lesion size
xxx<-exp(rnorm(N.pat,52.17,34.377))
tl_b.mean <-52.17
tl_b.sd <- 34.377
min.tl_b<- 7
max.tl_b<- 156 
xxx<-exp(rnorm(N.pat,tl_b.mean,tl_b.sd))
xxx <- pmax(pmin(xxx, max.tl_b), min.tl_b)
tl_b_wins<-rep(xxx,rep(n.times,N.pat))
#
#
# DIGARFST - time since initial diagnosis
DIGARFST.mean <-11.3756
DIGARFST.sd <- 25.86883
min.DIGARFST<- 0.1314
max.DIGARFST<- 189.83 
xxx<-exp(rnorm(N.pat,DIGARFST.mean,DIGARFST.sd))
xxx <- pmax(pmin(xxx, max.DIGARFST), min.DIGARFST)
DIGARFST<-rep(xxx,rep(n.times,N.pat))
#
#
#
#
#
# Create a data frame 'data.base' by combining all baseline variables
data.base<-data.frame(id=id,time=time,plannedtrt=plannedtrt,AGE=AGE,AGEGR1N=AGEGR1N,
                      ICNSMTFL= df$ICNSMTFL,PCHEMFL=df$PCHEMFL,MICNSFL=df$MICNSFL,ALKFDA=df$ALKFDA,
                      DIAGSTE=df$DIAGSTE,ECOG=ECOG,ECOGGR1=ECOGGR1,ECOGGR2=ECOGGR2,
                      RACEGR1=df$RACEGR1,SMKHIS=SMKHIS, SMKHISGR1=SMKHISGR1,
                      HISCLSE=HISCLSE,HISCLSEGR1=HISCLSEGR1,HISCLSEGR2=HISCLSEGR2,
                      SEX=df$SEX,
                      DIAGINITGR1=DIAGINITGR1,DIAGINITGR2=DIAGINITGR2,DIAGINITGR3=DIAGINITGR3,
                      STRATRAN=STRATRAN,
                      tl_b_wins=tl_b_wins,PRACTDS=df$PRACTDS,PRRADYN=df$PRRADYN,PRRBRYN=df$PRRBRYN,
                      LIASE=LIASE,DIGARFST=DIGARFST,REGION1N=REGION1N)
#
#
# Quadratic and cubic time
data.base$time2 <- data.base$time^2
#data.base$time2[data.base$time == -1] <- 0
data.base$time3 <- data.base$time^3
#data.base$time3[data.base$time == -1] <- 0
# Quadratic Age
data.base$AGE2 <- data.base$AGE^2
#
#
str(data.base)
table(table(data.base$id))
summary(data.base$AGE)
summary(data.base$tl_b_wins)
table(data.base$plannedtrt[1:1000],data.base$id[1:1000])
#
#
#
#
#
# Create time-dependent variables
# update the time-dependent variables
#
#
# Define a predictive function for binary outcomes
predict.binary.rfc<-function(beta=c(-1,1.4,2.1,-1.5),x.values=c(1,0.5,0.8,1))
{
  # Calculate the probability of success using a logistic model
  pp<-exp(sum(beta*x.values))
  pp<-pp/(1+pp)
  # Generate a binary outcome (0 or 1) based on the calculated probability
  return(rbinom(1,1,pp))
}
# Call the predict.binary.rfc function to generate a binary outcome
predict.binary.rfc()
#
#
# Define a predictive function for continuous outcomes
predict.continous.rfc<-function(beta=c(-1,1.4,2.1,-1.5),x.values=c(1,0.5,0.8,1),sigma=sqrt(21.75))
{
  # Calculate the mean of the continuous outcome using a linear model
  mm<-sum(beta*x.values)
  # Generate a continuous outcome based on the calculated mean and specified standard deviation (sigma)
    return(abs(rnorm(1,mm,sigma)))
}
# Call the predict.continuous.rfc function to generate a continuous outcome
predict.continous.rfc()
#
#
# Time dependent data
# Create a data frame 'data.td' to store time-dependent variables
#
# progtd: progression status Y/N (1/0)
# ecogtdGR2: ECOG 01 
# icprogtd: intracranial progression status Y/N (1/0)
# tltd_wins: target lesion size 
# trttd: treatment over time 01
# deathtd: death Y/N (1/0)
# censoringtd_new: censoring Y/N (1/0)
#
data.td<-data.frame(progtd=NA,ecogtdGR2=NA,icprogtd=NA,tltd_wins=NA,trttd=NA, deathtd=NA,censoringtd_new=NA)
# 
#
# Create time-dependent data for each unique patient ID
for (ind in unique(id))
# for (ind in 1:2)
  { 
#  ind<-20
  # Extract the data for the current patient ID (ind)
  xx<-data.base[data.base$id==ind,]
  
  # Initialize empty vectors for time-dependent variables
  ecogtdGR2<-vector("numeric",n.times)
  tltd_wins<-vector("numeric",n.times)
  progtd<-vector("numeric",n.times)
  deathtd<-vector("numeric",n.times) 
  icprogtd<-vector("numeric",n.times)     
  censoringtd_new<-vector("numeric",n.times)
  trttd<-vector("numeric",n.times)
  
  # Set the initial values of time-dependent variables as their baseline values, or as 0 at baseline
  ecogtdGR2[1]<-xx$ECOGGR1[1] 
  tltd_wins[1]<-xx$tl_b_wins[1]
  progtd[1]<-0
  deathtd[1]<-0
  icprogtd[1]<-0
  censoringtd_new[1]<-0
  trttd[1]<-xx$plannedtrt[1]

  # Update time-dependent variables for subsequent time points
  for (ttt in 2:n.times)
    {
    
    # Define coefficients and input values for progression prediction
    bb.prog.td<-c(-0.5,0.007,-0.001,-0.121,0.001,-0.715,0.923,0.423,0.271, 2.064, 0.009,-0.7)
    ww.prog.td<-c(1,xx$time[ttt],xx$time2[ttt],xx$AGE[ttt],xx$AGE2[ttt], xx$RACEGR1[ttt],xx$ECOGGR1[ttt],xx$SMKHISGR1[ttt],
                  xx$MICNSFL[ttt],icprogtd[ttt-1],tltd_wins[ttt-1],trttd[ttt-1])
    
    # Check if the lengths of beta and x values match
    if (length(bb.prog.td)!=length(ww.prog.td)) stop("Prog: Beta and x have different lengths")
    
    # Predict progression status for the current time point
    progtd[ttt]<-predict.binary.rfc(bb.prog.td,ww.prog.td)
    
    # Apply deterministic knowledge from initial trial
    if (progtd[ttt]==0 & xx$plannedtrt[ttt]==0) trttd[ttt]<-0
    # progtd is absorbing
    if (progtd[ttt-1]==1) progtd[ttt]<-1
    #
    #
    bb.ecogtdGR2<-c(-5,-0.003,-0.043,0.0003,0.398,1.130,-0.407,-0.126,-0.538,-0.254,0.839,-0.017,
                    0.088,-0.219,-0.399,0.458,0.155,0.834,-0.479,0.001,0.535,0.199,-0.502,-1)
    ww.ecogtdGR2<-c(1,xx$time[ttt],xx$AGE[ttt],xx$AGE2[ttt],xx$RACEGR1[ttt],xx$ECOGGR1[ttt],xx$SMKHISGR1[ttt],
                    xx$LIASE[ttt]==1,xx$LIASE[ttt]==2,xx$LIASE[ttt]==3,xx$DIAGINITGR3[ttt]==1,xx$DIAGINITGR3[ttt]==2,
                    xx$DIAGINITGR3[ttt]==3,xx$MICNSFL[ttt],xx$PRRADYN[ttt],
                    xx$STRATRAN[ttt]==1,xx$STRATRAN[ttt]==2,xx$STRATRAN[ttt]==3,icprogtd[ttt-1],tltd_wins[ttt-1],
                    progtd[ttt-1],progtd[ttt],trttd[ttt-1],ecogtdGR2[ttt-1])
    # Generate 3 categories of ECOG
    if (ecogtdGR2[ttt-1]==0)
    {
      ecogtdGR2[ttt]<-predict.binary.rfc(bb.ecogtdGR2,ww.ecogtdGR2)
    }
    else
    if (ecogtdGR2[ttt-1]==1)
    {
      ecogtdGR2[ttt]<-1+predict.binary.rfc(bb.ecogtdGR2,ww.ecogtdGR2)
    }
    else
      ecogtdGR2[ttt]<-2
   if (length(bb.ecogtdGR2)!=length(ww.ecogtdGR2)) stop("ECOG: Beta and x have different lengths")
    #
    #
    bb.icprog.td<-c(1.6,-0.020,-0.164,0.001,-0.866,0.452,-1.158,1.892,0.538,2.776,
                    -0.032,0.103,1.303,-0.276,-1.036,-17.041,-0.002,-3.523,3.589,-1.370)
    ww.icprog.td<-c(1,xx$time[ttt],xx$AGE[ttt],xx$AGE2[ttt],xx$RACEGR1[ttt],xx$MICNSFL[ttt],xx$PRRADYN[ttt],
                    xx$STRATRAN[ttt]==1,xx$STRATRAN[ttt]==2,xx$STRATRAN[ttt]==3,xx$LIASE[ttt]==1,xx$LIASE[ttt]==2,
                    xx$LIASE[ttt]==3,xx$DIAGINITGR3[ttt]==1,xx$DIAGINITGR3[ttt]==2,xx$DIAGINITGR3[ttt]==3,
                    tltd_wins[ttt-1],progtd[ttt-1],progtd[ttt],trttd[ttt-1])
    if (length(bb.icprog.td)!=length(ww.icprog.td)) stop("icprog: Beta and x have different lengths")
    icprogtd[ttt]<-predict.binary.rfc(bb.icprog.td,ww.icprog.td)
    # icprogtd is absorbing
    icprogtd[ttt]<-ifelse(icprogtd[ttt-1]==1,1,icprogtd[ttt])
    #
    #
    bb.tl.wins.td<-c(-15.324,0.741,-0.031,0.0004,0.013,-0.00004,-0.204,-0.119,0.102,0.156,0.110,0.336,
                     0.090,-0.685,-0.263,0.202,-0.144,0.186,-0.189,0.120,0.508,-2.740,-0.301,0.976,
                     0.455,0.716,3.424,-0.017,-1.226)
    ww.tl.wins.td<-c(1,xx$time[ttt],xx$time2[ttt],xx$time3[ttt],xx$AGE[ttt],xx$AGE2[ttt],xx$SEX[ttt],xx$RACEGR1[ttt],
                     xx$ECOGGR1[ttt], xx$LIASE[ttt]==1,xx$LIASE[ttt]==2,xx$LIASE[ttt]==3,xx$DIAGINITGR3[ttt]==1,
                     xx$DIAGINITGR3[ttt]==2,xx$DIAGINITGR3[ttt]==3,xx$SMKHISGR1[ttt],xx$MICNSFL[ttt],xx$PRRADYN[ttt],
                     xx$STRATRAN[ttt]==1,xx$STRATRAN[ttt]==2,xx$STRATRAN[ttt]==3,progtd[ttt-1],
                     ecogtdGR2[ttt-1],tltd_wins[ttt-1],icprogtd[ttt-1],ecogtdGR2[ttt],progtd[ttt],icprogtd[ttt],trttd[ttt-1])
    if (length(bb.tl.wins.td)!=length(ww.tl.wins.td)) stop("TL: Beta and x have different lengths")
    tltd_wins[ttt]<-predict.continous.rfc(bb.tl.wins.td,ww.tl.wins.td)
    #
    #
    bb.trt.td<-c(-9,-0.1,0.200,-0.001,0.560,-0.504,1.634,-0.330,0.017,-0.887,-0.517,
                 0.110,-0.005,2.704,-1.544)
    ww.trt.td<-c(1,xx$time[ttt],xx$AGE[ttt],xx$AGE2[ttt],xx$RACEGR1[ttt],xx$STRATRAN[ttt]==1,
                 xx$STRATRAN[ttt]==2,xx$STRATRAN[ttt]==3,tltd_wins[ttt-1],xx$DIAGINITGR3[ttt]==1,
                 xx$DIAGINITGR3[ttt]==2,xx$DIAGINITGR3[ttt]==3,tltd_wins[ttt],icprogtd[ttt-1],
                 icprogtd[ttt]) # add progtime , progtime2
    if (length(bb.trt.td)!=length(ww.trt.td)) stop("trt: Beta and x have different lengths")
    trttd[ttt]<-predict.binary.rfc(bb.trt.td,ww.trt.td)
    # no switching when planned trt =1
    if (xx$plannedtrt[ttt]==1) trttd[ttt]<-1
    # switching only possible at or after progression
    if (progtd[ttt] == 0) trttd[ttt] <- 0 
    # trttd is absorbing
    if (trttd[ttt-1]==1){
      trttd[ttt]<-1
    }
    
    #
    #
    bb.death.td<-c(-27,0.2,0.492,-0.003,0.593,-0.142,0.775,1.383,-0.114,0.019,-0.254,1.540,-1.232,
                   -1.372,-1.332,1.028,0.390,-1.154,-0.006,-0.003)
    ww.death.td<-c(1,xx$time[ttt],xx$AGE[ttt],xx$AGE2[ttt],xx$SMKHISGR1[ttt],xx$STRATRAN[ttt]==1,
                   xx$STRATRAN[ttt]==2,xx$STRATRAN[ttt]==3,ecogtdGR2[ttt-1],xx$LIASE[ttt]==1,xx$LIASE[ttt]==2,xx$LIASE[ttt]==3,xx$DIAGINITGR3[ttt]==1,
                   xx$DIAGINITGR3[ttt]==2,xx$DIAGINITGR3[ttt]==3,ecogtdGR2[ttt],trttd[ttt-1],trttd[ttt],
                   tltd_wins[ttt-1],tltd_wins[ttt])# progtime, progtime2
    if (length(bb.death.td)!=length(ww.death.td)) stop("Death: Beta and x have different lengths")
    deathtd[ttt]<-predict.binary.rfc(bb.death.td,ww.death.td)
    # deathtd is absorbing
    deathtd[ttt]<-ifelse(deathtd[ttt-1]==1,1,deathtd[ttt])
    # deathtd only possible at or after progression
    deathtd[ttt]<-ifelse(progtd[ttt]==0,0,deathtd[ttt])## AT
    #
    #
    bb.cens.td<-c(-18,0.45,-0.004,0.0001,-0.126,0.001,-0.230,0.243,-0.973,0.946,1.632,-1.407,
                  0.868,1.913,-0.802,-0.066,-1.620,14.593,-13.418,-0.003,0.016,-1.084,1.827,1.041,
                  -0.544,15.850,-13.486)
    ww.cens.td<-c(1,xx$time[ttt],xx$time2[ttt],xx$time3[ttt],xx$AGE[ttt],xx$AGE2[ttt],xx$LIASE[ttt]==1,
                  xx$LIASE[ttt]==2,xx$LIASE[ttt]==3, xx$SEX[ttt],xx$RACEGR1[ttt],
                  xx$ECOGGR1[ttt],xx$SMKHISGR1[ttt],xx$PRRADYN[ttt],xx$STRATRAN[ttt]==1,xx$STRATRAN[ttt]==2,xx$STRATRAN[ttt]==3,
                  icprogtd[ttt-1],icprogtd[ttt],tltd_wins[ttt-1],tltd_wins[ttt],progtd[ttt-1],progtd[ttt],
                  ecogtdGR2[ttt],ecogtdGR2[ttt-1],trttd[ttt-1],trttd[ttt])
    if (length(bb.cens.td)!=length(ww.cens.td)) stop("Cens: Beta and x have different lengths")
    censoringtd_new[ttt]<-predict.binary.rfc(bb.cens.td,ww.cens.td)
    # censoring is absorbing
    censoringtd_new[ttt]<-ifelse(censoringtd_new[ttt-1]==1,1,censoringtd_new[ttt])
    #
    #
  }
  # Create a data frame 'xx.td' containing time-dependent variables for the current patient
    xx.td<-cbind(ecogtdGR2,tltd_wins,progtd,deathtd,icprogtd,censoringtd_new,trttd)#,progtime
  data.td<-rbind(data.td,xx.td)
}

dim(data.td)
data.td<-data.td[-1,]  
#
#
#
#
#
# Merge Baseline df (data.base) and time-dependent df (data.td)
data.all<-cbind(data.base,data.td)
data.all$ecogtdGR2<- factor(data.all$ecogtdGR2)
dim(data.base)
#
#
#
str(data.all)
xx<-ls()
ls()
#
#
table(table(data.all$id))
#
# Data management
library(dplyr)
# Remove observations after death or censoring
data1 <- data.all %>%
  group_by(id) %>%
  arrange(id, time) %>% 
  filter(row_number() <= which.max(cumsum(deathtd == 1 | censoringtd_new == 1) == 1)) %>% # Filter rows up to the first event occurrence
  mutate(
    deathind = ifelse(sum(deathtd)== 1, 1, 0), # Create a death indicator variable
    censoringind = ifelse(sum(censoringtd_new) == 1, 1, 0))%>%
  
  mutate(
    deathtd = ifelse(censoringtd_new == 1 & deathind == 0, NA, deathtd), # Update 'deathtd' values to NA when 'censoringtd_new' is 1 and 'deathind' is 0
    censoringtd_new = ifelse(!is.na(censoringtd_new) & !is.na(deathtd) & censoringtd_new == 1 & deathtd == 1, 0, censoringtd_new))%>% # Update 'censoringtd_new' values to 0 when death is 1
  ungroup()
#
#
#
# Create xotd (switching status over time)
# Create xotime (time to switch)
# Create progtime (time to progression)
data2<- data1 %>%
  dplyr::group_by(id)%>%
  dplyr:: mutate(
    maxtime= max(time),
    xotd = ifelse(plannedtrt==1,0,
                 ifelse(plannedtrt==0 & trttd==0,0,
                        ifelse(plannedtrt==0 & trttd==1,1,NA))),
    first_switch_index = ifelse(sum(xotd == 1) > 0, min(which(xotd == 1)), NA),
    xotd = ifelse(!is.na(first_switch_index) & row_number() > first_switch_index, NA, xotd),
    xotime = ifelse(!is.na(first_switch_index), time[first_switch_index], maxtime),
    first_prog_index = ifelse(sum(progtd == 1) > 0, min(which(progtd == 1)), NA),
    progtime = ifelse(!is.na(first_prog_index), time[first_prog_index], maxtime))%>%
  ungroup()
#
#
# Create maxvisit_cens, deathoverall_new
data3<- data2 %>%
  dplyr::group_by(id)%>%
  dplyr::mutate(
    cens_visit=ifelse(is.na((which(xotd %in% 1)[1])),max(time),(which(xotd %in% 1)[1]-2)),
    cens_new= ifelse(cens_visit==time,1,ifelse(cens_visit<time,NA,0)),
    maxvisit_cens= ifelse(cens_visit==-1,0,cens_visit),
    deathoverall_new= ifelse(is.na(cens_new),NA,
                             ifelse(is.na((which(deathtd %in% 1)[1])),0,
                                    ifelse((which(deathtd %in%1)[1]-1)<= maxvisit_cens,1,0))))%>%
  dplyr::ungroup()
#
#
# Create binary time-varying ecogtdGR2 from ecogtdGR1
#
data3$ecogtdGR1 <- ifelse(data3$ecogtdGR2 %in% c(1, 2), 1, 0)
#
#
# Removing unnecessary variables
# data4 <- data3 %>%
#   select(-AGEGR1N,-ICNSMTFL,-PCHEMFL, -ALKFDA,-DIAGSTE,-ECOG,-ECOGGR2,-SMKHIS,
#          -HISCLSE,-HISCLSEGR1,-HISCLSEGR2,-DIAGINITGR1,-DIAGINITGR2,-tl_b_wins,
#          -PRACTDS,-PRRBRYN,-DIGARFST,-REGION1N,censoringind,-first_switch_index,-first_prog_index,
#          -cens_visit,-cens_new)
#
#
data4 <- data3[, !names(data3) %in% c(
  "AGEGR1N", "ICNSMTFL", "PCHEMFL", "ALKFDA", "DIAGSTE", "ECOG", "ECOGGR2", "SMKHIS",
  "HISCLSE", "HISCLSEGR1", "HISCLSEGR2", "DIAGINITGR1", "DIAGINITGR2", "tl_b_wins",
  "PRACTDS", "PRRBRYN", "DIGARFST", "REGION1N", "censoringind", "first_switch_index",
  "first_prog_index", "cens_visit", "cens_new"
)]

#
# Data Exploration:

# Here I do check for switchers: people with treatment 0 and with treatment 1.
# I look at a table which presents for each participant how many 0 and how many 1 are present
# A switcher needs to have "0" treatment and "1" treatment.
# If the result below is NULL that there are no switchers in the data set.
#
tt<-table(data4$trttd,data4$id)
str(tt)
tt[,1:200]
colnames(tt)[tt[1,]!=0 & tt[2,]!=0]
#
#
# how long certain patients stay in your data set
table(table(data4$id))
#
#
tt.death<-table(data4$deathtd,data4$id)
tt.death[,1:30]
tt.censoring<-table(data4$censoringtd_new,data4$id)
tt.censoring[,1:30]
table(data4$deathtd)
tt.progression<-table(data4$progtd,data4$id)
tt.progression[,1:30]
#
table(data4$plannedtrt,data4$trttd)
table(data4$plannedtrt,data4$deathtd)
table(data4$plannedtrt,data4$censoringtd_new)
table(data4$plannedtrt,data4$progtd)
table(data4$plannedtrt,data4$icprogtd)
#
#
#
table(data4$deathtd)
table(data4$censoringtd_new)
table(data4$censoringtd_new, data4$deathtd, useNA = "always")
table(data4$progtd)
table(data4$icprogtd)
#
#
#
summary(data4$deathtd)
summary(data4$censoringtd_new)
summary(data4$progtd)
summary(data4$tltd_wins)
summary(data4$trttd)
#
#
# SAVE** 
df_ipcw_pp<- data4
save(df_ipcw_pp, file = "df_ipcw_pp.rda")
#
#
#
