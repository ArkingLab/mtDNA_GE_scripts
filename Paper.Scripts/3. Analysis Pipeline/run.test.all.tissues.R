
library(yangR)
library(analyzeR)
library(pbapply)

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
  stop("One argument must be supplied (fileID).", call.=FALSE)
 }

num<-args[1]
num <- as.numeric(num)

#########################
### get the tissue ######
#########################
load('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/R_objects/all_tissuesv8.rds')

char <- all_tissues[num]

# set the directory

setwd(paste0('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/', char))

# check if the file exists
if (!file.exists('full.frame.rds')){print(paste0('full.frame.rds does not exist for ', char))} else{
load('full.frame.rds')


# 1.13.2020 changing to inverse normal transformation, adding race
# 1.31.2020 changing cohort to a factor....
full.frame.ct.corr$COHORT <- as.factor(full.frame.ct.corr$COHORT)

#########################
### do the lm ###########
#########################
num.svs <- length(grep('SV', colnames(full.frame.ct.corr)))

# right now you're using all SVs discovered.
# covnames <- c(paste0("SV", 1:sva.results$n.sv), paste0("Genotyping.PC", 1:3), "SEX", "COHORT")
if(length(levels(droplevels(as.factor(full.frame.ct.corr$COHORT)))) == 1){
	covnames <- c(paste0("SV", 1:num.svs), "SEX", "PC1", "PC2", "PC3", "RACE", "smtsisch") 
} else if(length(levels(droplevels(as.factor(full.frame.ct.corr$SEX)))) == 1){
	covnames <- c(paste0("SV", 1:num.svs), "COHORT", "PC1", "PC2", "PC3", "RACE", 'smtsisch') 
} else{covnames <- c(paste0("SV", 1:num.svs), "SEX", "COHORT", "PC1", "PC2", "PC3", "RACE", 'smtsisch') }

print(covnames)

full.frame.ct.corr$mtDNA.inv.norm <- inv.norm.transform(full.frame.ct.corr$mtDNA_adjust_AGE)
my.list <- get.vars(full.frame.ct.corr, permutevar = 'mtDNA.inv.norm', covnames = covnames)

# grep('ENSG00000213866.3', colnames(full.frame.ct.corr))
# expr <- full.frame.ct.corr[,grep('ENSG00000213866.3', colnames(full.frame.ct.corr))]

run_lm_default2 <- function(expr, cov, SCORE, omit.outlier = T, outlier_sd = 3) {
  expr <- as.numeric(expr)
  # expr <- scale(expr) # YOU SCALED EXPRESSION HERE! don't scale.  (but would that really make a diff?)

  # DO AN INVERSE NORMAL TRANSFORMATION!!!
  expr <- inv.norm.transform(expr)
  expr_cov <- cbind(SCORE, expr, cov)
  
  
  # # Run lm() normal procedure
  lm.fit <- lm(expr ~ ., data = expr_cov)
  lm.fit.summary <- summary(lm.fit)
  
  # Get corr
  cor_expr_score <-
    cor(expr, SCORE)
  
  # Capture p-val, etc.
  lm.res_ <-
    as.data.frame(t(coef(lm.fit.summary)['SCORE',]))
  intercept <- coef(lm.fit)[1]
  names(intercept) <- NULL
  lm.res_ <-
    cbind(data.frame(intercept),
          lm.res_,
          t(confint(lm.fit)['SCORE', ]),
          cor_expr_score)
  
  return(as.matrix(lm.res_))
}

run_lm2 <- function(expr, cov, SCORE, omit.outlier = T, method='default') {
  res <- switch (method,
                 "default" = run_lm_default2(expr, cov, SCORE, omit.outlier),
                 "two-stage" = run_lm_two_stage(expr, cov, SCORE, omit.outlier)
  )
  return(res)
}

run.all.lms2 <- function(tx_expr, cov, gene.ids, SCORE, omit.outlier = T, num.cores = 10)
{
  require(pbapply)
  lm.res <-
    pblapply(tx_expr,            # Expression vector list for `pbapply::pblapply`
             run_lm2,             # This function
             cov = cov,           # Covariate matrix, as desribed above
             SCORE = SCORE,       # PRS
             omit.outlier = omit.outlier,
             method = 'default',# Choose between 'default' or 'two-stage'
             cl = num.cores)      # Number of cores to parallelize over
  
  lm.res <- simplify2array(lm.res, higher=F)
  rownames(lm.res) <-
    c('intercept',
      'beta',
      'SE',
      't_value',
      'pval',
      'conf.low',
      'conf.high',
      'corr.rho')
  colnames(lm.res) <- gene.ids
 lm.res <- as.data.frame(t(lm.res))
  return(lm.res)
  # Sort results by p-value
  # No sorting! This will mess up the permutation saving
  # lm_res.sort <- lm.res[order(lm.res$pval), ]
  # return(lm_res.sort)
}

lm_res.sort <- run.all.lms2(my.list[[1]], my.list[[2]], my.list[[3]], my.list[[4]], omit.outlier = T, num.cores = 10)


# my.list <- get.vars(full.frame.ct.corr, permutevar = 'mtDNA_adjust_AGE', covnames = covnames)
# lm_res.sort2 <- run.all.lms(my.list[[1]], my.list[[2]], my.list[[3]], my.list[[4]], omit.outlier = T, num.cores = 10)


####### create a place for images ########
# if(dir.exists('images') == F) {dir.create('images')}

png(file = paste0('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/pval.qqplot.images/', char, '.pval.qqplot.png'), width = 700, height = 400)
pval_qqplot(lm_res.sort$pval, title = paste0('QQ-plot for ', char)) 
dev.off() 

# load gene key and merge
load('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/R_objects/v8.gene_key.rds')
lm_res.sort$gene_id <- rownames(lm_res.sort)
with.gene <- merge(lm_res.sort, gene_key, by = 'gene_id')
with.gene <- with.gene[order(with.gene$pval),]

to.look <- dplyr::select(with.gene, symbol, beta, pval)
print(head(to.look))

save(with.gene, file = 'with.gene.rds')
}
