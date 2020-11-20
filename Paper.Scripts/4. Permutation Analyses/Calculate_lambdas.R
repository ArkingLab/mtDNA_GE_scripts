# 2.20.2020
# Reran 4.20.2020
# Reran 5.13.2020
# Changed to 1000 permutations 11.14.2020

# The purpose of this script is to attempt to find an appropriate genomic inflation factor cutoff to determine which tissues have significant global association with mtDNA-CN
# This will also get the permutation cutoff for significance for each tissue!

# Lambda calculations:
# setwd
setwd('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8')

# load libraries
# install stuff
library(devtools)
source('/dcs01/arking/arkinglab/users/syang/libLoad.R')

# load full list of tissues
load('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/R_objects/all_tissuesv8.rds')

# load gene key
load('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/R_objects/v8.gene_key.rds')



all.lambdas <- as.data.frame(matrix(nrow = 1, ncol = 2))
colnames(all.lambdas) <- c('Lambda', 'Tissue')

for (i in 1:length(all_tissues)){
        char <- all_tissues[i]
        print(paste0("On tissue ", char))
        wdir <- paste0('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/new.permutation.method/', char)
        if(file.exists(paste0('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/new.permutation.method/', char, '/perm.1'))==F)
        {
                print('This directory does not exist')} else{
                setwd(wdir)
                for(x in 1:1000){
                        filename <- paste0('perm.', x)
					      pvals <- as.data.frame(fread(filename))
					      pvals=pvals[,2]
					      # print(paste0('On permutation ', x))
                        chisq <- qchisq(1-pvals,1)
                        lambda<-median(chisq,na.rm=TRUE) / 0.456
                        all.lambdas <- rbind(all.lambdas, c(lambda, char))
                }
                # png(paste0('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/new.permutation.method/images/', char, '.png'))
                # pval_qqplot.permutes(perm.frame[,1], perm.frame, char)
                # dev.off()
        }
}


# Kidney should NOT be in here!!! You ran it cause somebody asked you to....
all.lambdas = na.omit(all.lambdas)
all.lambdas = subset(all.lambdas, Tissue != 'Kidney - Cortex')

all.lambdas <- all.lambdas[order(all.lambdas$Lambda, decreasing = T), ]

all.lambdas$Lambda <- as.numeric(all.lambdas$Lambda)

all.lambdas$Lambda[nrow(all.lambdas) * 0.05]

save(all.lambdas, file = '/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/R_objects/all.lambdas.rds')

load('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/R_objects/lambda.frame.rds')

sig.tiss <- subset(lambda.frame, Lambda > all.lambdas$Lambda[nrow(all.lambdas) * 0.05])
nonsig.tiss <- subset(lambda.frame, Lambda < all.lambdas$Lambda[nrow(all.lambdas) * 0.05])
# 48 tissues were tested, including blood


cutoffs <- as.data.frame(matrix(nrow = 1, ncol = 2))
colnames(cutoffs) <- c('Tissue', 'Cutoff')

for (i in 1:length(all_tissues)){
  char <- all_tissues[i]
  wdir <- paste0('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/new.permutation.method/', char)
  if(file.exists(paste0('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/new.permutation.method/', char, '/perm.1'))==F)
  { print('This directory does not exist')} else{
    pvals = numeric()
    setwd(wdir)
    for(x in 1:1000){
      filename <- paste0('perm.', x)
      cutoff <- as.data.frame(fread(filename))
      pvals=c(pvals, min(cutoff[,2]))
    }
    pvals = pvals[order(pvals)]
    pvals[5]
    cutoffs = rbind(cutoffs, c(char, pvals[5]))
    print(paste('Finished with', char))
  }
}

cutoffs = na.omit(cutoffs)
perm.cutoffs = cutoffs
save(perm.cutoffs, file = '/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/new.permutation.method/cutoffs.rds')





