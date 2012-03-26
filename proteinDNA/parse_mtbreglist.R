# Script to parse the MtbRegList flat file obtained from 
# Pierre-E. Jacques [jacquespe@gis.a-star.edu.sg] on March 20, 2012 at
# reistett@ohsu.edu

#OHSU

#Laptop
setwd('~/schoolDB/Dropbox/thesis_work/')

#Read in flat file
dat.raw <- read.table('./data/protein-DNA/MtbRegList-litt_R1-1.txt',
                      head=F, stringsAsFactors=F, quote='', sep='\t')

#Not interested in TSS, just regulators and promoters. Select those.
dat.regs <- dat.raw[dat.raw[,1]%in%c("Regulateur", "Promoteur"),c(2,3)]

#Don't care about consensus sequence, just want reg-target pairs. So
#filter out those.
consensus.idx <- which(toupper(dat.regs[,2])=="CONSENSUS")
dat.filtered <- dat.regs[-consensus.idx,]



#Convert gene names (e.g. devR) to Rv IDs using biomaRt
library(biomaRt)
bac <- useMart('bacteria_mart_13', dataset='myc_30_gene')

#Get all the names in the first column
gene.names.1.idx <- grepl("rv", dat.filtered[,1], ignore.case=TRUE)
gene.names.1 <- dat.filtered[!gene.names.1.idx, 1]

gene.ids.1 <- getBM(attributes=c("tuberculist", "external_gene_id"), 
                    filters="external_gene_id", 
                    values=unique(gene.names.1), mart=bac)

#Not all the names returned a value, check which ones
unique(gene.names.1)[which(!(toupper(unique(gene.names.1)) %in% toupper(gene.ids.1[,2])))]
#[1] "FurA"    "Unknown" "SigH" 

#Manually look up FurA and SigH in TBDB
#SigH: Rv3223c
#FurA: Rv1909c

#Add them to IDs
gene.ids.1 <- rbind(gene.ids.1, c("Rv3223c", "sigH"))
gene.ids.1 <- rbind(gene.ids.1, c("Rv1909c", "furA"))

#Use the annotations to update all the IDs in dat.filtered to Rv IDs
for (i in 1:nrow(gene.ids.1)){
  id.idx <- which(toupper(dat.filtered[,1]) == toupper(gene.ids.1[i,2]))
  dat.filtered[id.idx,1] <- gene.ids.1[i,1]
}

#Get all the names in the first column
gene.names.2.idx <- grepl("rv", dat.filtered[,2], ignore.case=TRUE)
gene.names.2 <- dat.filtered[!gene.names.2.idx, 2]

gene.ids.2 <- getBM(attributes=c("tuberculist", "external_gene_id"), 
                    filters="external_gene_id", 
                    values=unique(gene.names.2), mart=bac)

#Not all the names returned a value, check which ones
unique(gene.names.2)[which(!(toupper(unique(gene.names.2)) %in% toupper(gene.ids.2[,2])))]
#[1] "fura"  "bfra"  "rrs"   "siga"  "sigb"  "sigH"  "sigh"  "trxB1" "trxb2"

#Manually look up in TBDB
#rrs not a gene in tbdb
gene.ids.2 <- rbind(gene.ids.2, c("Rv1909c", "furA"))
gene.ids.2 <- rbind(gene.ids.2, c("Rv1876", "bfrA"))
gene.ids.2 <- rbind(gene.ids.2, c("Rv2703", "sigA"))
gene.ids.2 <- rbind(gene.ids.2, c("Rv2710", "sigB"))
gene.ids.2 <- rbind(gene.ids.2, c("Rv3223c", "sigH"))
gene.ids.2 <- rbind(gene.ids.2, c("Rv1471", "trxB1"))
gene.ids.2 <- rbind(gene.ids.2, c("Rv3913", "trxB2"))

for (i in 1:nrow(gene.ids.2)){
  id.idx <- which(toupper(dat.filtered[,2]) == toupper(gene.ids.2[i,2]))
  dat.filtered[id.idx,2] <- gene.ids.2[i,1]
}
                                  
#Filter our pairs that didn't get updated to Rv IDs
good.1.idx <- grep("rv", dat.filtered[,1], ignore.case=TRUE)
dat.1clean <- dat.filtered[good.1.idx,]

good.2.idx <- grep("rv", dat.1clean[,2], ignore.case=TRUE)

mtbreglist.pairs <- dat.1clean[good.2.idx,]

colnames(mtbreglist.pairs) <- c("regulator", "target")


