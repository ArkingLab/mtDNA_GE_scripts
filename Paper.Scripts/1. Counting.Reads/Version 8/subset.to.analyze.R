library(data.table)

phenotypes <- as.data.frame(fread('/dcs01/arking/arkinglab/resources/GTeX/dbGaP_GTEx_phs000424.v8.p2/files/GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt'))

blood.rnaseq <- subset(phenotypes, SMTSD == 'Whole Blood' & SMAFRZE == 'RNASEQ')
blood.wgs <- subset(phenotypes, SMTSD == 'Whole Blood' & SMAFRZE == 'WGS') # no run info here....

library(yangR)

# filtering of people
blood.wgs$subjid <- make.subjids(blood.wgs$SAMPID)
blood.rnaseq$subjid <- make.subjids(blood.rnaseq$SAMPID)
both <- merge(blood.wgs, blood.rnaseq, by = 'subjid')
both$SMNABTCHD2 <- as.POSIXct(strptime(both$SMNABTCHD.x,format="%m/%d/%Y"))
after.2013.blood <- subset(both, SMNABTCHD2>'2013-01-01 EST')
write.csv(after.2013.blood$SAMPID, file = '/dcs01/arking/arkinglab/resources/GTeX/full.cramlist.lst')