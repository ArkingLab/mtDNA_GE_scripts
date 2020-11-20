# 5.11.2020 
# 11.15.2020 --> rerun with 1000 

# This script will get correlations between blood estimates and tissue estimates for all genes that were significant in blood
# It will also get correlations for random genes, to see what the null is.

# Load blood info
load('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/Whole Blood/with.gene.rds')
blood.with.gene = with.gene

# get pval cutoff
load('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/new.permutation.method/cutoffs.rds')
# pval.blood = as.numeric(perm.cutoffs$Cutoff[which(perm.cutoffs$Tissue == 'Whole Blood')])
pval.blood = 2.38e-6
sig.blood = subset(blood.with.gene, pval < pval.blood)

# omit negatives
# sig.blood.noneg = subset(sig.blood, beta > 0)

# load all tissues
load('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/R_objects/all_tissuesv8.rds')

# Create a dataframe for the results
results = as.data.frame(matrix(nrow = 1, ncol = 3))
colnames(results) = c('Tissue', 'Num.sig.blood', 'Correlation')

for(i in 1:54)
{
	char <- all_tissues[i]
	print(char)
	setwd(paste0('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/', char))

	if (!file.exists('with.gene.rds')){print(paste0('with.gene.rds does not exist for ', char))} else{
		load('with.gene.rds')
		both.sig = merge(sig.blood, with.gene, by = 'gene_id')

#		ggplot(both.sig, aes(beta.x, beta.y)) + xlab('Beta in Blood') + ylab(paste0('Beta in ', char)) + geom_point()

		cor.add = cor(both.sig$beta.x, both.sig$beta.y, method = 'spearman')
		to.add = c(char, nrow(both.sig), cor.add)

		results = rbind(results, to.add)
	}
}

spear.results = na.omit(results)
save(spear.results, file = '/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/R_objects/spear.results.rds')

# make a dataframe with 100 columns and 721 rows, each row is a randomly selected gene from whole blood
for.permute.genes = as.data.frame(matrix(nrow = nrow(sig.blood), ncol = 100))

set.seed(1)
for(i in 1:1000){
	rand.genes = sample(blood.with.gene$gene_id, nrow(sig.blood))
	for.permute.genes[,i] = rand.genes
}

perm.results = as.data.frame(matrix(nrow = 1, ncol = 3))
colnames(perm.results) = c('Tissue', 'Spear.Correlation', 'Permutation')

for(i in 1:54)
{
	char <- all_tissues[i]
	# print(char)
	setwd(paste0('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/', char))

	if (!file.exists('with.gene.rds')){print(paste0('with.gene.rds does not exist for ', char))} else{
		load('with.gene.rds')
		print(paste0('On tissue ', char))
		for (x in 1:1000){
			sig.blood = subset(blood.with.gene, gene_id %in% for.permute.genes[,x])
			both.sig = merge(sig.blood, with.gene, by = 'gene_id')
			cor.add = cor(both.sig$beta.x, both.sig$beta.y, method = 'spearman')
			to.add = c(char, cor.add, x)
			perm.results = rbind(perm.results, to.add)
		}
	}
}


# this can be WAY more efficient
# for (x in 1:1000){
# 	print(paste0('On permutation ', x))
# 	sig.blood = subset(blood.with.gene, gene_id %in% for.permute.genes[,x])
# 	for(i in 1:54)
# 	{
# 		char <- all_tissues[i]
# 		# print(char)
# 		setwd(paste0('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/', char))

# 		if (!file.exists('with.gene.rds')){print(paste0('with.gene.rds does not exist for ', char))} else{
# 			load('with.gene.rds')
# 			both.sig = merge(sig.blood, with.gene, by = 'gene_id')
# 			cor.add = cor(both.sig$beta.x, both.sig$beta.y, method = 'spearman')
# 			to.add = c(char, cor.add, x)

# 			perm.results = rbind(perm.results, to.add)
# 		}
# 	}

# }

perm.results = na.omit(perm.results)
perm.beta.results = perm.results
save(perm.beta.results, file = '/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/R_objects/perm.beta.results.rds')

# # load spearman results
# load('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/R_objects/spear.results.rds')
# spear.results$Spear.Correlation = as.numeric(spear.results$Correlation)
# perm.beta.results$Spear.Correlation = as.numeric(perm.beta.results$Spear.Correlation)

# library(ggplot2)
# # subset out char, see if actual result is diff from null
# length(unique(na.omit(perm.beta.results$Tissue)))

# unique(perm.beta.results$Tissue)

# list.chars = character()

# for(i in 1:length(unique(na.omit(perm.beta.results$Tissue))))
# {
# 	char <- unique(na.omit(perm.beta.results$Tissue))[i]
# 	tiss.perms = subset(perm.beta.results, Tissue == char)
# 	tiss.result = subset(spear.results, Tissue == char)
# 	# p = ggplot(tiss.perms, aes(Spear.Correlation)) + geom_density() + geom_rug() 
# 	# p + geom_text(data = tiss.result, aes(Spear.Correlation, y = 3, label = 'Observed correlation', vjust = -0.5), angle = 90, col = 'black') + geom_vline(xintercept = tiss.result$Spear.Correlation, lty = 2, colour = "gray50") + xlab('Distribution of null spearman correlations') + ggtitle(char)

# 	# How to define significance? 
# 	# If abs(tiss.result$Spear.Correlation) > abs(max(tiss.perms$Spear.Correlation))?
# 	top = tiss.perms$Spear.Correlation[order(tiss.perms$Spear.Correlation)]
	
# 	small = min(top)
# 	big = max(top)
# 	# small = top[length(top) * 0.025]
# 	# big = top[length(top) - (length(top) * 0.025)]

# 	if(tiss.result$Spear.Correlation > big | tiss.result$Spear.Correlation < small)
# 	{
# 		list.chars = c(list.chars, char)
# 	}
# }



# list.chars = character()

# for(i in 1:length(unique(na.omit(perm.beta.results$Tissue))))
# {
# 	char <- unique(na.omit(perm.beta.results$Tissue))[i]
# 	tiss.perms = subset(perm.beta.results, Tissue == char)
# 	tiss.result = subset(spear.results, Tissue == char)
# 	p = ggplot(tiss.perms, aes(Spear.Correlation)) + geom_density() + geom_rug() 
# 	p + geom_text(data = tiss.result, aes(Spear.Correlation, y = 3, label = 'Observed correlation', vjust = -0.5), angle = 90, col = 'black') + geom_vline(xintercept = tiss.result$Spear.Correlation, lty = 2, colour = "gray50") + xlab('Distribution of null spearman correlations') + ggtitle(char)

# 	# How to define significance? 
# 	# If abs(tiss.result$Spear.Correlation) > abs(max(tiss.perms$Spear.Correlation))?
# 	top = abs(tiss.perms$Spear.Correlation[order(tiss.perms$Spear.Correlation)])
# 	top = top[order(top)]
# 	top[950]
# 	# small = min(top)
# 	# big = max(top)
# 	# small = top[length(top) * 0.025]
# 	# big = top[length(top) - (length(top) * 0.025)]

# 	if(abs(tiss.result$Spear.Correlation) > top[950])
# 	{
# 		list.chars = c(list.chars, char)
# 	}
# }


# list chars

# save(list.chars, file = '/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/R_objects/sig.spear.tiss.rds')
# load('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/R_objects/perm.beta.results.rds')
# load('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/R_objects/spear.results.rds')
# spear.results$Spear.Correlation = as.numeric(spear.results$Correlation)
# perm.beta.results$Spear.Correlation = as.numeric(perm.beta.results$Spear.Correlation)

# spear.results$Set = 'Observed'
# perm.beta.results$Set = 'Permuted'

# all.results = dplyr::bind_rows(spear.results, perm.beta.results)

# ggplot(all.results, aes(Spear.Correlation, fill = Set)) + geom_rug() + geom_density()+ geom_text(data = spear.results, aes(Spear.Correlation, y = 2, label = Tissue, vjust = -0.5), angle = 90, col = 'black') + geom_vline(xintercept = spear.results$Spear.Correlation, lty = 2, colour = "gray50") 

# sig = subset(spear.results, Tissue %in% list.chars)
# library(htmlTable)
# sig = dplyr::select(sig, Tissue, Spear.Correlation)

# htmlTable(sig,
#           header =  c('Tissue', 'Spearman correlation'),
#           rnames = rep('', nrow(sig)))
        



