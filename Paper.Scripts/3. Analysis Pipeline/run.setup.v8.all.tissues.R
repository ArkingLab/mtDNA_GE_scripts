args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
  stop("One argument must be supplied (fileID).", call.=FALSE)
 }

num<-args[1]
num <- as.numeric(num)

setwd('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8')

# load libraries
# install stuff
library(devtools)
# install_github('syyang93/yangR')
# install_github('syyang93/analyzeR')
library(pbapply)
library(data.table)
library(yangR)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("edgeR")
# BiocManager::install("sva")
library(edgeR)
library(sva)

source('/dcs01/arking/arkinglab/users/syang/libLoad.R')

#############################
####### start analysis ######
#############################

# get tissue list
load('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/R_objects/all_tissuesv8.rds')

char <- all_tissues[num]

load(paste0('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/counts.rds.files/', char, '.rds')) # something is wrong with adipose specific counts

# ### Create directory for the tissue you are working on
if(dir.exists(char) == F) {dir.create(char)}

# set the directory to save in.

setwd(paste0('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/', char))

# get mtDNA-CN from blood info!
load('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/R_objects/final.mito.rds')

cn.only <- dplyr::select(final.mito, SUBJID, mtDNA, mtDNA_adjust_AGE, COHORT, PC1, PC2, PC3, SEX, AGE, RACE, TRISCHD)

tissue.counts <- specific.counts
rownames(tissue.counts) <- tissue.counts$Name
tissue.counts <- tissue.counts[,-c(1,2)]

tissue.counts_t <- as.data.frame(t(tissue.counts))
tissue.counts_t$SUBJID <- make.subjids(rownames(tissue.counts_t))

which(tissue.counts_t$SUBJID %in% cn.only$SUBJID) %>% length 

cn.only <- cn.only[which(cn.only$SUBJID %in% tissue.counts_t$SUBJID),]
phenotypes <- as.data.frame(fread('/dcs01/arking/arkinglab/resources/GTeX/dbGaP_GTEx_phs000424.v8.p2/files/GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt'))
head(phenotypes)
colnames(phenotypes) <- tolower(colnames(phenotypes))

#####################################################
########### change this line for tissue!! ###########
#####################################################
tiss.pheno <- subset(phenotypes, smtsd == char & smafrze == 'RNASEQ')

tiss.pheno$SUBJID <- make.subjids(tiss.pheno$sampid)
full.frame <- merge(tiss.pheno, cn.only, by = 'SUBJID')

w.only.txt <- tissue.counts_t[which(tissue.counts_t$SUBJID %in% cn.only$SUBJID),]
nrow(w.only.txt)

# if fewer than 50, then don't bother!!
if(nrow(w.only.txt) <= 50){print('Fewer than 50 samples, not run!')} else{

#############################################################################################
# Genes were selected based on expression thresholds of >0.1 TPM in at least 20% of samples #
#############################################################################################
tissue.counts_t <- w.only.txt
no_runs <- tissue.counts_t
no_runs$SUBJID <- NULL

no_runs <- apply(no_runs, 2, as.numeric)
rownames(no_runs) <- rownames(tissue.counts_t)

print("filtering transcripts")

zeroes <- colSums(no_runs <= 0.1) # greater than 0.1 TPMs
cutoff <- nrow(no_runs) * 0.8 # must have greater than 0.1 TPMs for this number of people

if (length(zeroes[which(zeroes > cutoff)]) != 0) {
	no_runs2 <- no_runs[, -which(zeroes > cutoff)]
	print('Number of genes removed:')
    print(ncol(no_runs) - ncol(no_runs2))
    print('Number of genes remaining:')
    print(ncol(no_runs2))
} else {
	no_runs2 <- no_runs
    print('No genes removed because of expression filter, genes remaining:')
    print(ncol(no_runs2))
} 

#############################
#### TMM normalization ######
#############################

print("Currently TMM normalizing count data")

library(edgeR)

save.no.runs <- no_runs2
no_runs2 <- as.data.frame(t(no_runs2))

y <- DGEList(counts=no_runs2)

# this shows you the normalization factors for each sample
normfactors <- calcNormFactors(y)
normfactors$sample

# this does the actual calculation!
tmm.normalized <- cpm(y)

look(tmm.normalized)
look(no_runs2)

tmm.normalized <- as.data.frame(t(tmm.normalized))
tmm.normalized$submitted_subject_id <- make.subjids(rownames(tmm.normalized))

save(tmm.normalized, file = 'tmm.normalized.rds')

# put the celltypes/phenotypes together! RNAseq covariates have already been added
full.frame.ct.corr <- merge(full.frame, tmm.normalized, by.x = 'SUBJID', by.y = 'submitted_subject_id')
full.frame <- full.frame.ct.corr

##########################################
########## get SVs #######################
##########################################
library(sva)

for_test <- full.frame[, grep("ENSG", colnames(full.frame))]
t.for.test <- t(for_test)

# full and null model specifications
full.mod <- model.matrix(~mtDNA_adjust_AGE, data = full.frame) # include all covariates you would like to protect from SVs in the full model
null.mod <- model.matrix(~1, data = full.frame) # have nothing in the null model

# t.for.test is my gene expression matrix
sva.results <- svaseq(t.for.test, full.mod, null.mod) 

svs <- as.data.frame(sva.results$sv)
colnames(svs) <- paste0('SV', 1:sva.results$n.sv)
# modsv <- cbind(lm_full, svs)
rownames(svs) <- full.frame$SUBJID
svs$SUBJID <- rownames(svs)

# outliers for SVs?
full.frame2 <- merge(full.frame, svs, by = 'SUBJID')

# there are outliers!
find_outliers <- function(svs)
{
    outlier <- c()
    all_outliers <- c()
    for (i in 1:(ncol(svs)-1)) {
        out <- find.outliers(svs[,i], sd = 4)
        outlier <- rownames(svs)[out]
        print(paste("outliers in SV", i, ":", sep = ""))
        print(outlier)
        all_outliers = c(all_outliers, outlier)
    }
    
    to.remove <- all_outliers
    to.remove <- unique(to.remove)
    print(paste0('total outliers: ', length(to.remove)))
    return(to.remove)
}

to.remove <- find_outliers(svs)

full.frame.save <- full.frame

##########################################
########## second iteration ##############
##########################################
if(length(to.remove) != 0){
removed.outliers <- full.frame[-which(full.frame$SUBJID %in% to.remove),]} else{removed.outliers <- full.frame}

nrow(removed.outliers)
full.frame <- removed.outliers
for_test <- full.frame[, grep("ENSG", colnames(full.frame))]
t.for.test <- t(for_test)

# full and null model specifications
full.mod <- model.matrix(~mtDNA_adjust_AGE, data = full.frame) # include all covariates you would like to protect from SVs in the full model
null.mod <- model.matrix(~1, data = full.frame) # have nothing in the null model
sva.results <- svaseq(t.for.test, full.mod, null.mod) 
svs <- as.data.frame(sva.results$sv)
colnames(svs) <- paste0('SV', 1:sva.results$n.sv)

rownames(svs) <- full.frame$SUBJID
svs$SUBJID <- rownames(svs)

to.remove <- find_outliers(svs)

full.frame.save <- full.frame

##########################################
########## third iteration ##############
##########################################
if(length(to.remove) != 0){
removed.outliers <- full.frame[-which(full.frame$SUBJID %in% to.remove),]} else{removed.outliers <- full.frame}

nrow(removed.outliers)
full.frame <- removed.outliers
for_test <- full.frame[, grep("ENSG", colnames(full.frame))]
t.for.test <- t(for_test)

# full and null model specifications
full.mod <- model.matrix(~mtDNA_adjust_AGE, data = full.frame) # include all covariates you would like to protect from SVs in the full model
null.mod <- model.matrix(~1, data = full.frame) # have nothing in the null model
sva.results <- svaseq(t.for.test, full.mod, null.mod)
svs <- as.data.frame(sva.results$sv)

colnames(svs) <- paste0('SV', 1:sva.results$n.sv)
rownames(svs) <- full.frame$SUBJID
svs$SUBJID <- rownames(svs)

# just to see how many outliers still left...
to.remove <- find_outliers(svs)
print(paste0(length(to.remove), ' outliers left after 3 iterations'))

save(svs, file = 'SVs.rds')

full.frame2 <- merge(full.frame, svs, by = 'SUBJID')
full.frame.ct.corr <- full.frame2
save(full.frame.ct.corr, file = 'full.frame.rds')
}
