# 11.11.2020
# Adapted from Dan Arking
# UKB analyses for neurodegenerative disease
# Make a table for manuscript

library(tidyverse)
library(lattice)
library(splines)
library(survival)
library(knitr)

## KEGG Pathways include Parkinsons and Alzheimers
## check for association in UKB
## Generate combined parkinsons and dementia (Neurodisease)
### check cases in each race group

# data2<-readRDS("/dcl01/arking/data/static/UKBiobank/GWAS_data/analyses/mtDNA_CN/ukbPheno_10072020_genotypes.rds")
data <-readRDS("/dcl01/arking/data/static/UKBiobank/GWAS_data/pheno/ukbPheno_09212020.rds")

# cool they are the same
# identical(data$arrayCN_PCAodd_m1, data2$arrayCN_PCAodd_m1)

# Dan already recoded race
##include ethnicity
# data$race<-fct_collapse(data$ethnicity,
#                         White = c("White","British","Irish","Any other white background"),
#                         Mixed = c("Mixed", "White and Black African", "White and Asian","Any other mixed background","White and Black Caribbean"),
#                         Asian = c("Asian or Asian British","Indian","Pakistani","Bangladeshi","Any other Asian background"),
#                         Black = c("Black or Black British","Caribbean","African","Any other Black background"),
#                         Chinese =c("Chinese"),
#                         Other = c("Other ethnic group"),
#                         Missing = c("Prefer not to answer","Do not know")
# )


# combine dementia and parkinsons --> no need to include Alzh because all individuals with Alzh have dementia
data <- data %>% mutate(Neurodisease=ifelse((Parkinson==1 | Dementia==1),1,0), incNeurodisease=ifelse((incParkinson==1 | incDementia==1),1,0), prvNeurodisease=ifelse((prvParkinson==1 | prvDementia==1),1,0),ttNeurodisease=pmin(ttDementia,ttParkinson))

# remove as incident if either disease is prevalent
data<- data %>% mutate(incNeurodisease=ifelse(prvNeurodisease==1,NA,incNeurodisease),ttNeurodisease=ifelse(prvNeurodisease==1,NA,ttNeurodisease))

t0<-data %>% filter(!(is.na(arrayCN_PCAodd_m1))) %>% group_by(Neurodisease,race) %>% summarise(Count=n())

t0    

# Need to use whites because the array metric was predicted using WES (mostly performed in whites), measure is not as reliable for other races
clean.data<-data %>% filter(race =="White" & used.in.pca.calculation==1)

########################################################################
## Functions for analysis 
## Focus on European ancestry, since all other groups <10 cases
## Only use unrelated individuals (used.in.pca.calculation==1)
## arrayCN_PCAodd_m1 is corrected for ns(age,df=2) + sex
## arrayCN_PCAodd_m2 is corrected for ns(age,df=2) + sex + Lymph + Platelet + Mono + NucRBC + Neutrophil + Baso + Eos
### cannot use too many covariates due to small number of cases
### Model 1 = age + sex + arrayCN_PCAodd_m1
### Model 1ex = age + sex + arrayCN_PCAodd_m1 (exluding cell type count outliers)
### Model 2 = age + sex + arrayCN_PCAodd_m2 (exluding cell type count outliers)
### Model 3 = model 2 + SmokingStatus
### Model 4 = model 3 + center
### Model 5 = model 4 + year + month

log_model<-function(pheno,dataset) { 
  pheno.model<-paste0("prv",pheno)
  model1<-c("age","sex","arrayCN_PCAodd_m1")
  
  lm1<-glm(reformulate(model1,pheno.model),data=dataset,family="binomial")
  lm1ex<-glm(reformulate(model1,pheno.model),data=subset(dataset,!(is.na(arrayCN_PCAodd_m2))),family="binomial")
  
  df = subset(dataset,!(is.na(arrayCN_PCAodd_m2)))
  control = length(which(df[, pheno.model] == 0))
  cases = length(which(df[, pheno.model] == 1))

  print(summary(lm1ex))
  
  #results
  m1<-c(format(coef(summary(lm1))["arrayCN_PCAodd_m1",],digits=3))
  m1ex<-c(format(coef(summary(lm1ex))["arrayCN_PCAodd_m1",],digits=3))
  
  HR = exp(coef(summary(lm1ex))["arrayCN_PCAodd_m1",1])

  lower = exp(coef(summary(lm1ex))["arrayCN_PCAodd_m1",1] - 1.96*(coef(summary(lm1ex))["arrayCN_PCAodd_m1",2]))
  upper = exp(coef(summary(lm1ex))["arrayCN_PCAodd_m1",1] + 1.96*(coef(summary(lm1ex))["arrayCN_PCAodd_m1",2]))
  
  pval = coef(summary(lm1ex))["arrayCN_PCAodd_m1",4]
  out = c(paste0('Prevalent ', pheno), HR, lower, upper, paste0(cases, '/', control), pval)
  return(out)
}

cox_model1<-function(pheno,dataset) { 
  pheno.model<-paste0("Surv(tt",pheno,",inc",pheno,")")
  model1<-c("age","sex","arrayCN_PCAodd_m1")
  model2<-c("age","sex","arrayCN_PCAodd_m2")
  model3<-c(model2,"SmokingStatus")
  model4<-c(model3,"Center")
  #model5<-c(model4,"year","month")
  lm1<-coxph(reformulate(model1,pheno.model),data=dataset)
  lm1ex<-coxph(reformulate(model1,pheno.model),data=subset(dataset,!(is.na(arrayCN_PCAodd_m2))))
  print(summary(lm1ex))
  
  #results
  m1<-c(format(coef(summary(lm1))["arrayCN_PCAodd_m1",],digits=3))
  m1ex<-c(format(coef(summary(lm1ex))["arrayCN_PCAodd_m1",],digits=3))
    
  n = lm1ex$n
  cases = lm1ex$nevent
  control = n - cases
  HR = summary(lm1ex)$conf.int["arrayCN_PCAodd_m1",1]
  lower = summary(lm1ex)$conf.int["arrayCN_PCAodd_m1",3]
  upper = summary(lm1ex)$conf.int["arrayCN_PCAodd_m1",4]
  pval = coef(summary(lm1ex))["arrayCN_PCAodd_m1",5]
  out = c(paste0('Incident ', pheno), HR, lower, upper, paste0(cases, '/', control), pval)
  return(out)
}

cox_model2<-function(pheno,dataset) { 
  pheno.model<-paste0("Surv(tt",pheno,",inc",pheno,")")
  model1<-c("age","sex","mtDNA_CN1")

  lm1<-coxph(reformulate(model1,pheno.model),data=dataset)
  lm1ex<-coxph(reformulate(model1,pheno.model),data=subset(dataset,!(is.na(mtDNA_CN2))))
  print(summary(lm1ex))
  
  #results
  m1<-c(format(coef(summary(lm1))["mtDNA_CN1",],digits=3))
  m1ex<-c(format(coef(summary(lm1ex))["mtDNA_CN1",],digits=3))
    
  n = lm1ex$n
  cases = lm1ex$nevent
  control = n - cases
  HR = summary(lm1ex)$conf.int["mtDNA_CN1",1]
  lower = summary(lm1ex)$conf.int["mtDNA_CN1",3]
  upper = summary(lm1ex)$conf.int["mtDNA_CN1",4]
  pval = coef(summary(lm1ex))["mtDNA_CN1",5]
  out = c(paste0('Incident ', pheno), HR, lower, upper, paste0(cases, '/', control), pval)
  return(out)
}

########################################################################

## Parkinsons
t1<-clean.data %>% filter(!(is.na(arrayCN_PCAodd_m1))) %>% group_by(Parkinson, prvParkinson) %>% summarise(Count=n())
t1

pd_prev = log_model("Parkinson",clean.data)
pd = cox_model("Parkinson",clean.data)

## Alzheimers
t2<-clean.data %>% filter(!(is.na(arrayCN_PCAodd_m1))) %>% group_by(Alzh, prvAlzh) %>% summarise(Count=n())
t2

ad = cox_model("Alzh",clean.data)
ad_prev = log_model('Alzh', clean.data)

## Dementia (excludes Alzh)
t3<-clean.data %>% filter(!(is.na(arrayCN_PCAodd_m1))) %>% filter(Alzh == 0) %>% group_by(Dementia, prvDementia) %>% summarise(Count=n())
t3

dem = cox_model("Dementia",clean.data %>% filter(Alzh==0))
dem_prev = log_model("Dementia",clean.data %>% filter(Alzh==0))

## Neurodisease
t5<-clean.data %>% filter(!(is.na(arrayCN_PCAodd_m1))) %>% group_by(Neurodisease,prvNeurodisease) %>% summarise(Count=n())
t5

nd = cox_model("Neurodisease",clean.data)
nd_prev = log_model("Neurodisease",clean.data)

show_inc = rbind(pd, ad, dem, nd)
show_prev = rbind(pd_prev, ad_prev, dem_prev, nd_prev)
colnames(show_inc) = c('Disease', 'Hazard ratio', 'Lower', 'Upper', 'Number of cases/controls', 'P-value')
colnames(show_prev) = c('Disease', 'Hazard ratio', 'Lower', 'Upper', 'Number of cases/controls', 'P-value')

save(show_inc, file = '/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/R_objects/UKB/show_inc.rds')
save(show_prev, file = '/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/R_objects/UKB/show_prev.rds')


# Time to event:
clean.data$ttAlz %>% na.omit %>% median
clean.data$ttParkinson %>% na.omit %>% median
clean.data$ttDementia %>% na.omit %>% median
clean.data$ttNeurodisease %>% na.omit %>% median

########################################################################
# Estimates for overlapping samples ####################################
########################################################################
full = clean.data[-which(is.na(clean.data$mtDNA_CN2)),]

give_z<-function(pheno,dataset) { 
  pheno.model<-paste0("Surv(tt",pheno,",inc",pheno,")")
  model1<-c("age","sex","arrayCN_PCAodd_m1")

  lm1ex<-coxph(reformulate(model1,pheno.model),data=subset(dataset,!(is.na(arrayCN_PCAodd_m2))))
  array_z<-c(format(coef(summary(lm1ex))["arrayCN_PCAodd_m1",4],digits=3))

  model1<-c("age","sex","mtDNA_CN1")
  
  lm1ex<-coxph(reformulate(model1,pheno.model),data=subset(dataset,!(is.na(mtDNA_CN2))))
  wes_z<-c(format(coef(summary(lm1ex))["mtDNA_CN1",4],digits=3))
 
  out = c(array_z, wes_z)
  return(out)
}

p_inc = give_z("Parkinson", full)
dem_inc = give_z("Dementia", full %>% filter(Alzh == 0))
alzh_inc = give_z("Alzh", full)
comb_inc = give_z("Neurodisease", full)

df = as.data.frame(rbind(p_inc, dem_inc, alzh_inc, comb_inc), stringsAsFactors = F)
colnames(df) = c('Array_z', 'WES_z')




give_z_prev <-function(pheno,dataset) { 
  pheno.model<-paste0("prv",pheno)

  model1<-c("age","sex","arrayCN_PCAodd_m1")
  lm1ex<-glm(reformulate(model1,pheno.model),data=subset(dataset,!(is.na(arrayCN_PCAodd_m2))),family="binomial")
  array_z<-format(coef(summary(lm1ex))["arrayCN_PCAodd_m1",3], digits = 3)

  model1<-c("age","sex","mtDNA_CN1")
  lm1ex<-glm(reformulate(model1,pheno.model),data=subset(dataset,!(is.na(arrayCN_PCAodd_m2))),family="binomial")
  wes_z<-format(coef(summary(lm1ex))["mtDNA_CN1",3], digits = 3)
 
  out = c(array_z, wes_z)
  return(out)
}

p_prev = give_z_prev("Parkinson", full)
dem_prev = give_z_prev("Dementia", full %>% filter(Alzh == 0))

rownames(df)
df = rbind(df, p_prev, dem_prev)
rownames(df)[c(4,5)] = c('p_prev', 'dem_prev')

df$WES_z = as.numeric(df$WES_z)
df$Array_z = as.numeric(df$Array_z)


lm(arrayCN_PCAodd ~ age, data = full) %>% summary
lm(resid.mtDNA ~ age, data = full) %>% summary

lm(arrayCN_PCAodd ~ sex, data = full) %>% summary
lm(resid.mtDNA ~ sex, data = full) %>% summary


p_prev = give_z_prev("Parkinson", full)
dem_prev = give_z_prev("Dementia", full %>% filter(Alzh == 0))




park_new = cox_model1("Parkinson", full)
park_old = cox_model2("Parkinson", full)

dem_new = cox_model1("Dementia", full %>% filter(Alzh == 0))
dem_old = cox_model2("Dementia", full %>% filter(Alzh == 0))

alzh_new = cox_model1("Alzh", full)
alzh_old = cox_model2("Alzh", full)

nd_new = cox_model1("Neurodisease", full)
nd_old = cox_model2("Neurodisease", full)

show_inc_all = as.data.frame(rbind(park_new, park_old, dem_new, dem_old, alzh_new, alzh_old, nd_new, nd_old))
colnames(show_inc_all) = c('Disease', 'Hazard ratio', 'Lower', 'Upper', 'Number of cases/controls', 'P-value')

save(show_inc_all, file = '/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/R_objects/UKB/show_inc_all.rds')


