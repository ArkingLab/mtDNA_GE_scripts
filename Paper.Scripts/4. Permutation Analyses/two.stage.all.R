### uses permutation method from jeff leek!
# 5.13.2020 # realized i wasn't scaling

# 2.20.2020
# This should create residuals and estimates for two-stage permutations for all tissues.
# inverse-normal transforms gene expression and mtDNA-CN metrics. 

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
	stop("One argument must be supplied (fileID).", call.=FALSE)
}

# load libs
library(analyzeR)
library(yangR)
library(magrittr)

# get the tissue
num<-args[1]
num <- as.numeric(num)
load('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/R_objects/all_tissuesv8.rds')
char <- all_tissues[num]

setwd(paste0('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/', char))

# load dataframe

if (!file.exists('full.frame.rds')){print(paste0('full.frame.rds does not exist for ', char))} else{
	load('full.frame.rds')

# 1. fit the alternative hypothesis, get residuals. 
	num.svs <- length(grep('SV', colnames(full.frame.ct.corr)))
	if(length(levels(droplevels(as.factor(full.frame.ct.corr$COHORT)))) == 1){
		covnames <- c(paste0("SV", 1:num.svs), "SEX", "PC1", "PC2", "PC3", "RACE", "smtsisch") 
	} else if(length(levels(droplevels(as.factor(full.frame.ct.corr$SEX)))) == 1){
		covnames <- c(paste0("SV", 1:num.svs), "COHORT", "PC1", "PC2", "PC3", "RACE", 'smtsisch') 
	} else{covnames <- c(paste0("SV", 1:num.svs), "SEX", "COHORT", "PC1", "PC2", "PC3", "RACE", 'smtsisch') }

	full.frame.ct.corr$mtDNA.inv.norm <- inv.norm.transform(full.frame.ct.corr$mtDNA_adjust_AGE)
	my.list <- get.vars(full.frame.ct.corr, permutevar = 'mtDNA.inv.norm', covnames = covnames)

# you're just going to run a SICK for loop. 

	expr.matrix <- my.list[[1]]
	cov <- my.list[[2]]
	gene.ids <- my.list[[3]]
	SCORE <- my.list[[4]]
	omit.outlier = T
	outlier_sd = 3

	resid.frame <- as.data.frame(matrix(nrow = nrow(expr.matrix), ncol = ncol(expr.matrix)))
	colnames(resid.frame) <- gene.ids
	rownames(resid.frame) <- full.frame.ct.corr$SUBJID

### INT OF GENES HERE ####
	expr.matrix2 <- apply(expr.matrix, 2, inv.norm.transform)
	expr.matrix <- as.data.frame(expr.matrix2)

##################
#### STEP 1 ######
##################

# getting residuals for copy number for each gene for each person: covariates = SVs 1-18, sex, cohort, race, and ischemic time

# mtDNA-CN 

# DONT SCALE RESIDUALS.

# check to make sure SVs are run on UNSCALED data
# could always use z scores if you wanted to directly compare
# comparing across genes, you're losing information about how much a gene is expressed

	for(i in 1:ncol(expr.matrix))
	{
		if(i %% 5000 == 0){print(paste('On gene', i))}
		expr <- as.numeric(expr.matrix[,i])
		# expr <- scale(expr)
		expr_cov <- cbind(SCORE, expr, cov)
		if (omit.outlier == T) {
			m <- mean(expr)
			s <- sd(expr)
			outliers <- which(expr > m + outlier_sd * s | expr < m - outlier_sd * s)
			if (length(outliers) != 0) {
				expr_cov$expr[outliers] <- NA
			}
		}
		lm.fit <- lm(expr ~ ., data = expr_cov, na.action = na.exclude)
		lm.resids <- resid(lm.fit)
		resid.frame[,i] <- lm.resids
	}
	resids.frame <- resid.frame

	save(resids.frame, file = paste0('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/', char, '/resids.frame.rds'))

##################
#### STEP 2 ######
##################
	estimates.frame <- as.data.frame(matrix(nrow = nrow(expr.matrix), ncol = ncol(expr.matrix)))
	colnames(estimates.frame) <- gene.ids
	rownames(estimates.frame) <- full.frame.ct.corr$SUBJID

# get estimates for each gene expression value, using only covariates and NO mtDNA-CN

# DONT SCALE ESTIMATES

	for(i in 1:ncol(expr.matrix))
	{
		if(i %% 5000 == 0){print(paste('On gene', i))}
		expr <- as.numeric(expr.matrix[,i])
		# expr <- scale(expr)

	# took out SCORE here.
		expr_cov <- cbind(expr, cov)
		if (omit.outlier == T) {
			m <- mean(expr)
			s <- sd(expr)
			outliers <- which(expr > m + outlier_sd * s | expr < m - outlier_sd * s)
			if (length(outliers) != 0) {
				expr_cov$expr[outliers] <- NA
			}
		}
		lm.fit <- lm(expr ~ ., data = expr_cov, na.action = na.exclude)
		lm.predicted <- predict(lm.fit)
		estimates.frame[,i] <- lm.predicted
	}

	fitted.frame <- estimates.frame

	save(fitted.frame, file = paste0('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/', char, '/fitted.frame.rds'))

##################
#### STEP 3 ######
##################
	# make the permutations frame!!

	order.for.permutes <- as.data.frame(matrix(nrow = nrow(full.frame.ct.corr), ncol = 100))

	set.seed(1)
	for(i in 1:100)
	{
		order.for.permutes[,i] <- sample(full.frame.ct.corr$SUBJID)
	}
	save(order.for.permutes, file = paste0('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/', char, '/order.for.permutes.rds'))
} 






