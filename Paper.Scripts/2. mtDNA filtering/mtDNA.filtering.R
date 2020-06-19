#### This script will take all.read.info.txt, and filter/QC it.
source('/dcs01/arking/arkinglab/users/syang/libLoad.R')

mtdna <- as.data.frame(fread('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/data.files/all.read.info.txt'))
ggplot(mtdna, aes(unaligned.reads)) + geom_density() + geom_vline(xintercept = (5*10^7), col = 'red')

# remove outliers for unaligned %s
nrow(mtdna)
mtdna <- mtdna[-find.outliers(mtdna$unaligned.reads),]
nrow(mtdna)

# make sure these are all from whole blood!
phenotypes <- as.data.frame(fread('/dcs01/arking/arkinglab/resources/GTeX/dbGaP_GTEx_phs000424.v8.p2/files/GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt'))
wgs.only <- subset(phenotypes, SMAFRZE == 'WGS')

both <- merge(mtdna, phenotypes, by.x = 'submitted_subject_id', by.y = 'SAMPID')
table(both$SMTSD)
# they are all from whole blood. we're good. 

# remove outliers past a certain cutoff for unaligned reads
mtdna <- mtdna[-which(mtdna$unaligned.reads > 5*10^7),]

mtdna$aligned.reads <- mtdna$Total.reads - mtdna$unaligned.reads
mtdna$mtDNA <- mtdna$mito.reads/mtdna$aligned.reads

mt <- mtdna

# get subject-level phenotype data:
subj.pheno <- as.data.frame(fread('/dcs01/arking/arkinglab/resources/GTeX/dbGaP_GTEx_phs000424.v8.p2/files/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt'))

mt$SUBJID <- make.subjids(mt$submitted_subject_id)
mt2 <- merge(mt, subj.pheno, by = 'SUBJID')

save.image('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/RDatas/make.mt.frame.12.17.19.RData')

mt2$RACE <- as.factor(mt2$RACE)
mt2$SEX <- as.factor(mt2$SEX)

lm(mtDNA ~ AGE, data = mt2) %>% summary
lm(mtDNA ~ SEX, data = mt2) %>% summary
lm(mtDNA ~ RACE, data = mt2) %>% summary

# add cell type info.
load('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/xCell/blood.xcell.blood.only.rds')

celltypes.t <- as.data.frame(t(xcell.blood.only))
colnames(celltypes.t) <- make.names(colnames(celltypes.t))
celltypes.t$SUBJID <- make.subjids(rownames(celltypes.t))

mt3 <- merge(mt2, celltypes.t, by = 'SUBJID')

celltypes <- colnames(celltypes.t)

for.form <- paste(celltypes, collapse=' + ')
for.names <- paste("\"", celltypes, collapse = ', ')

v = paste0('"', celltypes, '",')
# use cat in this case to see what's "really" there
# print will show the quotes escaped with backslashes
cat(v)

# run through all cell types, look for significant associations?

for(i in 1:(length(celltypes)-1)){
	type <- celltypes[i]
	form <- as.formula(paste0('mtDNA~', type))
	model <- lm(form, data = mt3) 
	print(type)
	print(coef(summary(model))[2,4])
}
# Neutrophils, Monocytes, very significant

#### backwards elimination #####
for.regress <- dplyr::select(mt3, mtDNA, celltypes)

for.regress$SUBJID <- NULL
library(MASS)
# Fit the full model 
full.model <- lm(mtDNA ~., data = for.regress)
# Stepwise regression model
step.model <- stepAIC(full.model, direction = "both", 
                      trace = FALSE)
summary(step.model)
# stepwise regression --> several celltypes still included in model

## pairs plots: ##
pairs(~mtDNA + CD4..T.cells + CD8..Tem + Neutrophils + NKT + Eosinophils + HSC + Macrophages.M1 + Macrophages.M2 + Megakaryocytes + Monocytes + Tregs + Th1.cells + Th2.cells, data = for.regress, lower.panel=twolines, diag.panel=mydiag.panel, upper.panel=panel.cor, label.pos=0.5, main="Correlation matrix of mtDNA and cell types", labels = c('mtDNA', 'CD4 cells', 'CD8 cells', 'Neutrophils', 'NKT', 'Eosinophils', 'HSC', 'Macrophages.M1', 'Macrophages.M2', 'Megakaryocytes', 'Monocytes', 'Tregs', 'Th1.cells', 'Th2.cells'))

lm(mtDNA ~ CD4..T.cells + CD8..Tem + Neutrophils + NKT + Eosinophils + HSC + Macrophages.M1 + Macrophages.M2 + Megakaryocytes + Monocytes + Tregs + Th1.cells + Th2.cells, data = for.regress) %>% summary

# remove CD8, Monocytes, and Th1cells (too much collinearity (>0.8 correlation))
lm(mtDNA ~ CD4..T.cells + Neutrophils + NKT + Eosinophils + HSC + Macrophages.M1 + Macrophages.M2 + Megakaryocytes + Tregs + Th2.cells, data = for.regress) %>% summary 

# remove nonsignificant cell types
lm(mtDNA ~ Neutrophils + NKT + Eosinophils + HSC + Megakaryocytes + Tregs, data = for.regress) %>% summary 

####################################
#### filter cell type outliers #####
####################################

find.outliers(mt3$CD4..T.cells)
find.outliers(mt3$Neutrophils)
find.outliers(mt3$NKT)
find.outliers(mt3$Eosinophils)
find.outliers(mt3$Macrophages.M1)
find.outliers(mt3$Macrophages.M2)
find.outliers(mt3$Megakaryocytes)
find.outliers(mt3$Tregs)
find.outliers(mt3$Th2.cells)

all.outliers <- c(find.outliers(mt3$CD4..T.cells), find.outliers(mt3$Neutrophils), find.outliers(mt3$NKT), find.outliers(mt3$Eosinophils), find.outliers(mt3$Macrophages.M1), find.outliers(mt3$Macrophages.M2), find.outliers(mt3$Megakaryocytes), find.outliers(mt3$Tregs), find.outliers(mt3$Th2.cells))

all.outliers <- unique(all.outliers)

# by eye: HSC and CD4..T.cells have extreme outliers.  remove these!
find.outliers(mt3$HSC)
find.outliers(mt3$CD4..T.cells)
mt4 <- mt3[-c(find.outliers(mt3$HSC), find.outliers(mt3$CD4..T.cells)),]

mt4$mtDNA_adjust_celltypes <- scale(resid(lm(mtDNA ~ CD4..T.cells + Neutrophils + NKT + Eosinophils + HSC + Megakaryocytes + Tregs, data = mt4)))

ggplot(mt4, aes(mtDNA, mtDNA_adjust_celltypes)) + geom_point() + xlab('No adjusting for celltypes') + ylab('Adjusting for celltypes')

library(splines)
lm(mtDNA_adjust_celltypes~ns(AGE,df=2), data = mt4) %>% summary
lm(mtDNA_adjust_celltypes~SEX, data = mt4) %>% summary

library(gridExtra)
grid.arrange(age, sex, ncol = 2)

ggplot(mt4, aes(RACE, mtDNA_adjust_celltypes)) + geom_boxplot() + stat_summary(fun.data = give.n, geom = 'text') + ylab('mtDNA ~ Cell types') 
ggplot(mt4, aes(COHORT, mtDNA_adjust_celltypes)) + geom_boxplot() + stat_summary(fun.data = give.n, geom = 'text') + ylab('mtDNA ~ Cell types') 

###############################
# drop surgery individual #####
###############################
nosurge <- subset(mt4, COHORT != 'Surgical')
mt4 <- nosurge

###############################
#### final mtdna-cn model #####
###############################

table(mt4$DTHHRDY) # i don't think you can use hardy score for categories of 14 people.....

lm(mtDNA ~ CD4..T.cells + Neutrophils + NKT + Eosinophils + HSC + Megakaryocytes + Tregs + ns(AGE,df=2) + SEX + COHORT + TRISCHD + as.factor(DTHHRDY), data = mt4) %>% summary

# take out HARDY
lm(mtDNA ~ CD4..T.cells + Neutrophils + NKT + Eosinophils + HSC + Megakaryocytes + Tregs + ns(AGE,df=2) + SEX + COHORT + TRISCHD, data = mt4) %>% summary

mt4$mtDNA_adjust_AGE <- scale(resid(lm(mtDNA ~ Neutrophils + HSC + Megakaryocytes + COHORT + TRISCHD + AGE + SEX, data = mt4)))
fit <- lm(mtDNA ~ Neutrophils + HSC + Megakaryocytes + COHORT + TRISCHD + AGE + SEX, data = mt4)
af <- anova(fit)
afss <- af$"Sum Sq"
print(cbind(af,PctExp=afss/sum(afss)*100))

# add genotyping PCs into this file!
genoPCs <- as.data.frame(fread('/dcs01/arking/arkinglab/resources/GTeX/dbGaP_GTEx_phs000424.v8.p2/files/GTEx_Analysis_2017-06-05_v8_WGS_VCF_files/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_support_files/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_20genotPCs.txt'))

mt5 <- merge(mt4, genoPCs, by.x = 'submitted_subject_id', by.y = 'FID')

# add race to model??
mt5$RACE <- factor(mt5$RACE, levels = c(3, 2, 1, 99))

# doesn't really make a difference
lm(mtDNA ~ CD4..T.cells + Neutrophils + NKT + Eosinophils + HSC + Megakaryocytes + Tregs + ns(AGE,df=2) + SEX + COHORT + TRISCHD + RACE, data = mt4) %>% summary

# PC1 and 2 don't really make a difference either
lm(mtDNA ~ CD4..T.cells + Neutrophils + NKT + Eosinophils + HSC + Megakaryocytes + Tregs + ns(AGE,df=2) + SEX + COHORT + TRISCHD + PC1 + PC2, data = mt5) %>% summary

lm(mtDNA ~ Neutrophils + HSC + Megakaryocytes + COHORT + TRISCHD + AGE + SEX) -> fit
af <- anova(fit)
afss <- af$"Sum Sq"
print(cbind(af,PctExp=afss/sum(afss)*100))


# write out corrected mtDNA-CN to an rds file
final.mito <- mt5
save(final.mito, file = '/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/R_objects/final.mito.rds')

######################
#### stopped here ####
######################
