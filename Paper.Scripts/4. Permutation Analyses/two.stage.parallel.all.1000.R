library(analyzeR)
library(yangR)
library(pbapply)

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
	stop("One argument must be supplied (fileID).", call.=FALSE)
}

num<-args[1]
num <- as.numeric(num)
load('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/R_objects/all_tissuesv8.rds')
char <- all_tissues[num]

if (!file.exists(paste0('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/', char, '/full.frame.rds'))){print(paste0('full.frame.rds does not exist for ', char))} else{
	
# figure out how to do permutations......
# 1. load the fitted frame
	load(paste0('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/', char, '/fitted.frame.rds'))

# 2. load the resids frame
	load(paste0('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/', char, '/resids.frame.rds')) 

# 3. load permutations frame
	load(paste0('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/', char, '/order.for.permutes1000.rds')) 

# 3.5. load the full frame 
	load(paste0('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/', char, '/full.frame.rds'))
# subset full frame on individuals in permutations frame

	full.frame.ct.corr <- full.frame.ct.corr[match(order.for.permutes[,1], full.frame.ct.corr$SUBJID),]
	full.frame.ct.corr <- full.frame.ct.corr[order(full.frame.ct.corr$SUBJID),]

	num.svs <- length(grep('SV', colnames(full.frame.ct.corr)))

# 4. get covariates, etc.
	num.svs <- length(grep('SV', colnames(full.frame.ct.corr)))
	if(length(levels(droplevels(as.factor(full.frame.ct.corr$COHORT)))) == 1){
		covnames <- c(paste0("SV", 1:num.svs), "SEX", "PC1", "PC2", "PC3", "RACE", "smtsisch") 
	} else if(length(levels(droplevels(as.factor(full.frame.ct.corr$SEX)))) == 1){
		covnames <- c(paste0("SV", 1:num.svs), "COHORT", "PC1", "PC2", "PC3", "RACE", 'smtsisch') 
	} else{covnames <- c(paste0("SV", 1:num.svs), "SEX", "COHORT", "PC1", "PC2", "PC3", "RACE", 'smtsisch') }

	full.frame.ct.corr$mtDNA.inv.norm <- inv.norm.transform(full.frame.ct.corr$mtDNA_adjust_AGE)
	
	# amygdala needs race taken out:
	 if(char == 'Brain - Amygdala')
	 {
	 covnames = covnames[-which(covnames == 'RACE')]
	 }

	my.list <- get.vars(full.frame.ct.corr, permutevar = 'mtDNA.inv.norm', covnames = covnames)

# 5. Make directory to save in
	if(dir.exists(paste0('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/new.permutation.method/', char)) == F){
		dir.create(paste0('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/new.permutation.method/', char))
	}
	setwd(paste0('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/new.permutation.method/', char))

# 6. Get permuted resids frame. --> set for a loop!

	for(i in 1:1000){
		samp.resids.frame <- resids.frame[match(rownames(resids.frame), order.for.permutes[,i]),]
		print(paste0('On permutation ', i))
# 7. Create null dataset
		null.data <- fitted.frame + samp.resids.frame        

		tx_expr <- null.data
		cov <- my.list[[2]]
		gene.ids <- my.list[[3]]
		SCORE <- my.list[[4]]

		lm.res <- run.all.lms(tx_expr, cov, gene.ids, SCORE, omit.outlier = T, num.cores = 10) # Amygdala cannot include RACE in the covnames...only three

		lm.res$pval
		print(min(lm.res$pval))
		write.table(lm.res$pval, file = paste0('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/new.permutation.method/', char, '/perm.', i))
		write.table(lm.res$beta, file = paste0('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/new.permutation.method/', char, '/perm.beta', i))
	}

}

