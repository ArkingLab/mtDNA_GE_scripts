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

library(EnvStats)
clean.data$resid.mtDNA_AGE = scale(resid(lm(resid.mtDNA ~ age, data = clean.data, na.action = na.exclude)))
clean.data$resid.mtDNA_SEX = scale(resid(lm(resid.mtDNA ~ sex, data = clean.data, na.action = na.exclude)))

# clean.data = clean.data[-which(is.na(clean.data$age)),]
clean.data = clean.data %>% filter(age>=40)
clean.data<-clean.data %>% mutate(agegroup=as.factor(cut(clean.data$age,seq(40,75,5),right=FALSE)))#,labels=c(1:7))))
age.plot<-clean.data %>%
  ggplot(aes(x=agegroup,y=resid.mtDNA)) +
  geom_boxplot() +
  stat_n_text() +
  labs(title="mtDNA-CN by age group (All)") + theme_classic() + ylab('Unadjusted mtDNA-CN') + xlab('Age group')
age.plot

clean.data$sex = factor(clean.data$sex, levels = c('Male', 'Female'))
sex.plot <- clean.data %>% ggplot(aes(x=as.factor(sex), y=resid.mtDNA)) + geom_boxplot() + labs(title="mtDNA-CN by sex") + theme_classic() + ylab('Unadjusted mtDNA-CN') + xlab('Biological sex')
sex.plot


lm(resid.mtDNA ~ sex, data = clean.data) %>% summary %>% coef
lm(resid.mtDNA ~ agegroup, data = clean.data) %>% summary %>% coef

library(patchwork)
age.plot + sex.plot

png(file = '/dcs01/arking/arkinglab/active/projects/GTeX/syang/ukb_ageplot.png')
age.plot + sex.plot
dev.off()

# Dan says add p-values

ggplot(clean.data, aes(age, resid.mtDNA)) + geom_point() + theme_classic() + geom_smooth(method = 'lm')