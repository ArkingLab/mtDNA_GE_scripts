###################################
####### load the RNAseq data ######
###################################

# load counts
counts <- as.data.frame(fread('/dcs01/arking/arkinglab/resources/GTeX/v8_GTEx_counts/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct'))

# you need the new phenotype file!!
phenotypes <- as.data.frame(fread('/dcs01/arking/arkinglab/resources/GTeX/dbGaP_GTEx_phs000424.v8.p2/files/GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt'))
head(phenotypes)

colnames(phenotypes) <- tolower(colnames(phenotypes))

all_tissues <- levels(as.factor(phenotypes$smtsd))

# save(all_tissues, file = '/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/R_objects/all_tissuesv8.rds)
load('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/R_objects/all_tissuesv8.rds')

char <- all_tissues[num]

save separate count files for each tissue.
look(counts)

for (i in 1:55)
{
    char <- all_tissues[i]
    print(char)
    tiss.counts <- subset(phenotypes, smtsd == char & smafrze == 'RNASEQ')
    genes.append <- counts[,1:2]
    specific.counts <- counts[,which(colnames(counts) %in% tiss.counts$sampid)]
    specific.counts <- cbind(genes.append, specific.counts)
    save(specific.counts, file = paste0('/dcs01/arking/arkinglab/active/projects/GTeX/syang/look.version8/counts.rds.files/', char, '.rds'))
}
