# 3.23.2020
# Modified 4.1.2020 --> Dan says to use T.test.pval instead of Rank.T.test.pval

# this script will perform permutation testing for GO enrichment (testing with blood first)
# It will also filter out pseudogenes before running the enrichments.

library(yangR)
library(analyzeR)
library(pbapply)
# if (!requireNamespace("BiocManager", quietly = TRUE))
  #  install.packages("BiocManager")

# BiocManager::install("qusage")

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
print(char)

# set the directory

setwd(paste0('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/', char))

# check if the file exists
if (!file.exists('with.gene.rds')){print(paste0('with.gene.rds does not exist for ', char))} else{
	
	load('with.gene.rds')
	load('/dcs01/arking/arkinglab/resources/ensembl/only.care.rds') # made in /dcs01/arking/arkinglab/resources/ensembl/try.read.R
	w.annot <- merge(only.care, with.gene, by = 'gene_id') # 23758 yay
	w.annot$symbol <- w.annot$symbol.y
	subset(w.annot, gene_type != ' ')
	# w.annot[which(w.annot$gene_type == ' Mt_rRNA'),]
	
	gene.types <- unique(w.annot$gene_type)
	pseudo.types <- gene.types[grep('pseudo', gene.types)]
	'%!in%' <- function(x,y)!('%in%'(x,y))
	with.gene <- subset(w.annot, gene_type %!in% pseudo.types)

	with.gene <- with.gene[order(with.gene$pval),]
	with.gene$abs.tscore <- abs(with.gene$t_value)

	# get ranked t-values
	with.gene$ranked.tvals <- order(with.gene$abs.tscore, decreasing = T) 

	# get gene sets
	library(qusage)
	library(data.table)

	# Kegg
	kegg.sets <- read.gmt('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/GO_enrich/MsigDB/c2.cp.kegg.v7.0.symbols.gmt')

	# GO
	go.sets <- read.gmt('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/GO_enrich/MsigDB/c5.all.v7.0.symbols.gmt')

	# TF
	tft.sets <- read.gmt('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/GO_enrich/MsigDB/c3.tft.v7.0.symbols.gmt')

	# function for testing
	perform.t.tests <- function(gene.sets, with.gene){
		signif.gene.sets <- as.data.frame(matrix(nrow = 1, ncol = 7))
		i <- 1
		colnames(signif.gene.sets) <- c('Gene.Set.Name', 'T.test.pval', 'Rank.t.test.pval', 'Num.genes.in.set', 'Beta', 'Confint.Upper', 'Confint.Lower')

		for(i in 1:length(gene.sets)){
			test.set <- gene.sets[[i]]
			set.name <- names(gene.sets[i])
			selected.indices <- which(with.gene$symbol %in% test.set)
			in.set <- with.gene[selected.indices,]
			out.set <- with.gene[-selected.indices,]

			# with.gene$Set = "Yes"
			# with.gene$Set[-selected.indices] = "No"			
			# library(ggplot2)
			# library(yangR)
			# ggplot(with.gene, aes(Set, abs.tscore)) + geom_boxplot() + stat_summary(fun.data = give.n, geom = 'text')
			if(nrow(in.set) < 5){ # if this is set to 2, for the ranked.t.stats some of the results will be 0!!!!
				# print(paste0("Not enough samples for ", set.name))
			} else{
				t.stats <- t.test(in.set$abs.tscore, out.set$abs.tscore)
				# if(t.stats$p.value < (0.05/length(gene.sets))){ # comment out if you want to see all gene sets....
					# print(set.name)
					# print(t.stats$p.value)
					ranked.t.stats <- t.test(in.set$ranked.tvals, out.set$ranked.tvals)
					# if(ranked.t.stats$p.value < 0.0005){ # get rid of this, its for troubleshooting
					# 	print(ranked.t.stats$p.value)
					# 	print(set.name)
					# }
					# print(ranked.t.stats$p.value)
					sig <- c(set.name, t.stats$p.value, ranked.t.stats$p.value, nrow(in.set), t.stats$estimate[1]-t.stats$estimate[2], t.stats$conf.int[1], t.stats$conf.int[2])
					signif.gene.sets <- rbind(signif.gene.sets, sig)
				# }
				i <- i + 1
				if(i %% 100 == 0){
				 	print(paste0('On term ', i , ' out of ', length(gene.sets)))
				}
			}
		}

		signif.gene.sets <- na.omit(signif.gene.sets)
		signif.gene.sets$T.test.pval <- as.numeric(signif.gene.sets$T.test.pval)

		signif.gene.sets <- signif.gene.sets[order(signif.gene.sets$T.test.pval, decreasing = F),]
		return(signif.gene.sets)
	}


	###############################
	#### KEGG, GO, TFT stuff ######
	###############################

	all.go.sets <- perform.t.tests(go.sets, with.gene)
	all.kegg.sets <- perform.t.tests(kegg.sets, with.gene)
	all.tft.sets <- perform.t.tests(tft.sets, with.gene)

	save(all.go.sets, file = 'all.go.rds')
	save(all.kegg.sets, file = 'all.kegg.rds')
	save(all.tft.sets, file = 'all.tft.rds')

	with.gene.permute <- with.gene
	set.seed(1)
 
	#### Get permutation cutoff for KEGG ####

	start <- Sys.time()
	all.min.kegg <- numeric()
	all.min.go <- numeric()
	all.min.tft <- numeric()

	for(p in 1:100){
		with.gene.permute$abs.tscore <- sample(with.gene.permute$abs.tscore)

		### kegg 
		permute.kegg <- perform.t.tests(kegg.sets, with.gene.permute)
		min.pval <- permute.kegg$T.test.pval[1]
		all.min.kegg <- c(all.min.kegg, min.pval)

		## go
		permute.go <- perform.t.tests(go.sets, with.gene.permute)
		min.pval <- permute.go$T.test.pval[1]
		all.min.go <- c(all.min.go, min.pval)

		### tft
		permute.tft <- perform.t.tests(tft.sets, with.gene.permute)
		min.pval <- permute.tft$T.test.pval[1]
		all.min.tft <- c(all.min.tft, min.pval)

		print(paste0('On permutation ', p))
	}
	end <- Sys.time()
	runtime <- end-start

	all.min.kegg <- all.min.kegg[order(all.min.kegg, decreasing = F)]

	print(paste0('Total runtime: ', runtime))

	### save permuted values
	all.perm.cutoff.forset <- data.frame(GO = all.min.go, KEGG = all.min.kegg, TFT = all.min.tft)
	
	all.perm.cutoff.forset$GO <- all.perm.cutoff.forset$GO[order(all.perm.cutoff.forset$GO, decreasing = F)]
	all.perm.cutoff.forset$KEGG <- all.perm.cutoff.forset$KEGG[order(all.perm.cutoff.forset$KEGG, decreasing = F)]
	all.perm.cutoff.forset$TFT <- all.perm.cutoff.forset$TFT[order(all.perm.cutoff.forset$TFT, decreasing = F)]

	save(all.perm.cutoff.forset, file = 'all.perm.cutoff.forset.rds')

	sig.kegg.results <- subset(all.kegg.sets, T.test.pval < all.perm.cutoff.forset$KEGG[nrow(all.perm.cutoff.forset) * 0.05])
	sig.go.results <- subset(all.go.sets, T.test.pval < all.perm.cutoff.forset$GO[nrow(all.perm.cutoff.forset) * 0.05])
	sig.tft.results <- subset(all.tft.sets, T.test.pval < all.perm.cutoff.forset$TFT[nrow(all.perm.cutoff.forset) * 0.05])	
	
	if(nrow(sig.kegg.results) != 0){sig.kegg.results$Set <- 'KEGG'}
	if(nrow(sig.go.results) != 0){sig.go.results$Set <- 'GO'}
	if(nrow(sig.tft.results) != 0){sig.tft.results$Set <- 'TFT'}		
	
	all.results <- plyr::rbind.fill(sig.kegg.results, sig.go.results, sig.tft.results)
	save(all.results, file = 'sig.pathways.nopseudo.permutecut.rds')
}

