# Re-updated on 4.7.2020
# This script will generate cell type compositions from Whole Blood data.
char <- 'Whole Blood'
setwd(paste0('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/', char))
load('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/R_objects/v8.gene_key.rds')
load(paste0('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/counts.rds.files/', char, '.rds'))


blood.counts <- specific.counts
rownames(blood.counts) <- blood.counts$Name
blood.counts <- blood.counts[,-c(1,2)]
blood.counts_t <- as.data.frame(t(blood.counts))

### make into format for xCell
rownames(blood.counts) <- gene_key$gene_id
blood.counts$symbol <- gene_key$symbol
# write.table(blood.counts, sep = '\t', file = '/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/xCell/blood.counts.txt', quote = F, row.names = gene_key$symbol)
# # this takes a while....

# try xCell in R:
devtools::install_github('dviraran/xCell')
library(xCell)
for.xcell <- blood.counts

####################################################
# filter out lowly expressed genes. ################
####################################################
no_runs <- blood.counts_t
no_runs$SUBJID <- NULL

no_runs <- apply(no_runs, 2, as.numeric)
rownames(no_runs) <- rownames(blood.counts_t)

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

b.counts <- as.data.frame(t(no_runs2))
b.counts$gene_id <- rownames(b.counts)
b.counts <- merge(b.counts, gene_key, by = 'gene_id')

####################################
# sum the genes that are the same ##
####################################
nrow(b.counts)
b.counts2 <- b.counts
b.counts2$gene_id <- NULL
aggregated <- aggregate(.~symbol, b.counts2, sum) # this takes a long time
nrow(aggregated) # looks like this worked! yay.

for.xcell <- aggregated
rownames(for.xcell) <- for.xcell$symbol
for.xcell$symbol <- NULL

xcell.all.cells <- xCellAnalysis(for.xcell) # only recognized 10,766 genes.

# # can run xCell using only a subset of cell types!
rownames(xcell.all.cells)
cell.types.use <- rownames(xcell.all.cells)[grep('CD', rownames(xcell.all.cells))]
cell.types.use <- c(cell.types.use, 'Neutrophils', 'Platelets', 'NKT', 'NK cells', 'Basophils', 'B-cells', 'Eosinophils', 'Erythrocytes', 'HSC', 'Macrophages', 'Macrophages M1', 'Macrophages M2', 'Mast cells', 'Megakaryocytes', 'Memory B-cells', 'Monocytes', 'Plasma cells', 'Tregs', 'Th1 cells', 'Th2 cells')

cell.types.use %!in% rownames(xcell.all.cells)

xcell.blood.only = xCellAnalysis(for.xcell, cell.types.use = cell.types.use) # 9298 genes used
 
save(xcell.all.cells, file = '/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/xCell/blood.xcell.all.cells.rds')
save(xcell.blood.only, file = '/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/xCell/blood.xcell.blood.only.rds')

