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


