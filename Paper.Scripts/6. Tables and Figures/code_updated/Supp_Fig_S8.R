# 11.8.2020

# This script will do PCA of brain tissues all together!
# load in all brain counts

load('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/R_objects/all_tissuesv8.rds')
library(analyzeR)
library(yangR)

all_counts = as.data.frame(fread('/dcs01/arking/arkinglab/resources/GTeX/v8_GTEx_counts/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct'))
phenotypes <- as.data.frame(fread('/dcs01/arking/arkinglab/resources/GTeX/dbGaP_GTEx_phs000424.v8.p2/files/GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt'))
brain = unique(phenotypes$SMTSD[grep('Brain', phenotypes$SMTSD)])

brain_only = subset(phenotypes, SMTSD %in% brain & SMAFRZE == 'RNASEQ')
brain_counts = all_counts[,which(colnames(all_counts) %in% brain_only$SAMPID)]
saveRDS(brain_counts, file = '/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/R_objects/brain_counts.rds')
brain_counts_t = as.data.frame(t(brain_counts))
colnames(brain_counts_t) = all_counts$Name
brain_counts_t$SAMPID = rownames(brain_counts_t)

tiss.only = dplyr::select(phenotypes, SAMPID, SMTSD)

w.tiss = merge(brain_counts_t, tiss.only, by = "SAMPID")

#############################################################################################
# Genes were selected based on expression thresholds of >0.1 TPM in at least 20% of samples #
#############################################################################################
tissue.counts_t <- w.tiss[,grep('ENSG', colnames(w.tiss))]

no_runs <- apply(tissue.counts_t, 2, as.numeric)
rownames(no_runs) <- rownames(brain_counts_t)

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

#############################
#### perform PCA ######
#############################
# maybe would be better to run a t-SNe (or a UMAP)

pca = prcomp(tmm.normalized, center = T, scale = F)
pcs = as.data.frame(pca$x)
pcs$tiss = w.tiss$SMTSD
saveRDS(pcs, file = '/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/R_objects/brain_pcs.rds')
ggplot(pcs, aes(PC1, PC2, col = tiss)) + geom_point() + theme_classic() + labs(col = "Tissue")
save(pcs, file = '/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/R_objects/brain_pcs_readable.rds')

# scale and center
pca = prcomp(tmm.normalized, center = T, scale = T)
pcs = as.data.frame(pca$x)
pcs$tiss = w.tiss$SMTSD
saveRDS(pcs, file = '/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/R_objects/brain_pcs_scale.rds')

tmm.normalized$submitted_subject_id <- make.subjids(rownames(tmm.normalized))
tmm.normalized$tissue = w.tiss$SMTSD
saveRDS(tmm.normalized, file = '/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/R_objects/brain_tmm.normalized.rds')

#############################
#### perform UMAP ######
#############################
library(umap)
labels = pcs$tiss
umap = umap(tmm.normalized)

plot = as.data.frame(umap$layout)
plot$tiss = labels

saveRDS(plot, file = '/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/R_objects/brain_umap.rds')

ggplot(plot, aes(V1, V2, col = tiss)) + geom_point()



