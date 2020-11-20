# 3.30.2020
# 4.12.2020
# 11.20.2020

# This script will look at GO enrichment across tissues (no pseudogenes) (using a permuted cutoff)

### This will look at GO enrichment across tissues 
# setwd
setwd('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8')

# load libraries
# install stuff
library(devtools)
source('/dcs01/arking/arkinglab/users/syang/libLoad.R')

# load full list of tissues
load('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/R_objects/all_tissuesv8.rds')

# Determine and save significant GO/KEGG for each tissue.

for (i in 1:length(all_tissues)){
	char <- all_tissues[i]
	print(paste0("Working on tissue, ", char))

	setwd(paste0('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/', char))
	# load all GO permutations
	if(file.exists('all.perm.cutoff.forset.rds')){
		load('all.perm.cutoff.forset.rds')
		all.perm = all.perm.cutoff.forset
		for(i in 1:12){
			load(paste0('all.perm.cutoff.forset_', (i*100), '.rds'))
			all.perm = rbind(all.perm, all.perm.cutoff.forset)	
		}

		# there are duplicates because 1-100 and 101-200 are the same seed. take them out!
		all.perm = all.perm[-which(duplicated(all.perm$GO)),]
		# you only need 1000 perms 
		all.perm = all.perm[1:1000,]
		# get tissue-specific cutoffs:
		go = all.perm$GO[order(all.perm$GO)]
		go.cutoff = go[nrow(all.perm)*0.05]

		kegg = all.perm$KEGG[order(all.perm$KEGG)]
		kegg.cutoff = kegg[nrow(all.perm)*0.05]

		tft = all.perm$TFT[order(all.perm$TFT)]
		tft.cutoff = tft[nrow(all.perm)*0.05]

		# load all results
		load('all.go.rds')
		sig.go = subset(all.go.sets, T.test.pval < go.cutoff)
		
		load('all.kegg.rds')
		sig.kegg = subset(all.kegg.sets, T.test.pval < kegg.cutoff)

		load('all.tft.rds')
		sig.tft = subset(all.tft.sets, T.test.pval < tft.cutoff)

		# make frame of sig results:
		all.results = rbind(sig.go, sig.kegg, sig.tft)
		if(nrow(all.results) != 0){
			all.results$Set = rep(c('GO', 'KEGG', 'TFT'), c(nrow(sig.go), nrow(sig.kegg), nrow(sig.tft)))
		}
		save(all.results, file = 'sig.pathways_1000.rds')
		}
}

# make a list of all significant genesets and terms
# save significant GO/KEGG for each tissue.

load('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/Whole Blood/sig.pathways_1000.rds')
big.sig.pathways <- all.results
big.sig.pathways$Tissue <- 'Whole Blood'

for (i in 1:(length(all_tissues)-1)){
	char <- all_tissues[i]
	print(paste0("Working on tissue, ", char))

	setwd(paste0('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/', char))
	# load all GO permutations
	if(file.exists('sig.pathways_1000.rds')){
		load('sig.pathways_1000.rds')
		if(nrow(all.results) != 0){
			all.results$Tissue = char
			big.sig.pathways = rbind(big.sig.pathways, all.results)
		}
	}
}

save(big.sig.pathways, file = '/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/R_code/Cross.tissue.look/big.sig.pathways.permutecut_1000.rds')

lookat <- as.data.frame(table(big.sig.pathways$Tissue))
lookat <- lookat[order(lookat$Freq, decreasing = T),]

lookgenes <- as.data.frame(table(big.sig.pathways$Gene.Set.Name))
lookgenes <- lookgenes[order(lookgenes$Freq, decreasing = T),]

load('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/R_code/Cross.tissue.look/big.sig.pathways.permutecut.rds')

# Most often-appearing term:
big.sig.pathways[grep('KEGG_RIBOSOME', big.sig.pathways$Gene.Set.Name),] # not in putamen

################################################################
## limit this so only looking at tissues with significant GIF ##
################################################################
load('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/R_objects/lambda.frame.rds')
load('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/R_objects/all.lambdas.rds')

lambda.cutoff <- all.lambdas$Lambda[nrow(all.lambdas) * 0.05] 

sig.tissues <- subset(lambda.frame, Lambda > lambda.cutoff)
nonsig.tissues <- subset(lambda.frame, Lambda < lambda.cutoff)

only.sig.tissues <- subset(big.sig.pathways, tissue %in% sig.tissues$Tissue)
lookgenes2 <- as.data.frame(table(only.sig.tissues$Gene.Set.Name))
lookgenes2 <- lookgenes2[order(lookgenes2$Freq, decreasing = T),]

##############################################
####### search full set for mito terms #######
##############################################
oxphos <- big.sig.pathways[grep('OXIDATIVE_PHOSPHORYLATION', big.sig.pathways$Gene.Set.Name),] # not in putamen
oxphos.tissues <- unique(oxphos$tissue)
length(oxphos.tissues) # 10

mito <- big.sig.pathways[grep('^GO_MITOCHONDRION$', big.sig.pathways$Gene.Set.Name),] # not in putamen
mito.tissues <- unique(mito$tissue)
length(mito.tissues) # 11

intersect(mito.tissues, oxphos.tissues) %>% length # there are only 5 overlaps :( )

union(mito.tissues, oxphos.tissues) %>% length # there are 16 shared

elk1 <- big.sig.pathways[grep('^SCGGAAGY_ELK1_02$', big.sig.pathways$Gene.Set.Name),] # not in putamen
elk1.tissues <- unique(elk1$tissue)
length(elk1.tissues)

union(mito.tissues, elk1.tissues) %>% length
union(oxphos.tissues, elk1.tissues) %>% length # there are 3 oxphos tissues not in elk1.tissues

which(mito.tissues %in% only.sig.tissues$tissue) %>% length # 11
which(oxphos.tissues %in% only.sig.tissues$tissue) %>% length # 10
which(elk1.tissues %in% only.sig.tissues$tissue) %>% length # 20


brain <- subset(big.sig.pathways, tissue == 'Brain - Putamen (basal ganglia)')


mito.tissues <- only.sig.tissues[grep("^GO_MITOCHONDRION$", only.sig.tissues$Gene.Set.Name),]
mito.tissues <- mito.tissues[order(mito.tissues$Rank.t.test.pval, decreasing = F),]





# make a list of all significant genes

load('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/Whole Blood/with.gene.rds')
load('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/Whole Blood/full.frame.rds')

big.with.gene.w.samps <- with.gene
big.with.gene.w.samps$Tissue <- 'Whole Blood'
big.with.gene.w.samps$N <- nrow(full.frame.ct.corr)

for (i in 1:(length(all_tissues)-1)){
	char <- all_tissues[i]
	print(paste0("Working on tissue, ", char))

	setwd(paste0('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/', char))
	# load all GO permutations
	if(file.exists('with.gene.rds')){
		load('with.gene.rds')
		if(nrow(with.gene) != 0){
			with.gene$Tissue = char
			load(paste0('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/', char, '/full.frame.rds'))
			with.gene$N <- nrow(full.frame.ct.corr)
			big.with.gene.w.samps = rbind(big.with.gene.w.samps, with.gene)
		}
	}
}

elk1_path = subset(big.sig.pathways, Gene.Set.Name == 'SCGGAAGY_ELK1_02')

elk1.only = subset(big.with.gene.w.samps, Tissue %in% elk1$Tissue)
elk1 = subset(elk1.only, symbol == 'ELK1')

save(elk1, file = '/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/R_objects/elk1.only.rds')

save(big.with.gene.w.samps, file = '/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/R_objects/big.with.gene.w.samps_1000.rds')



