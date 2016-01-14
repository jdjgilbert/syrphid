
### Import karyotype dataset

warning('NEED A DATA DICTIONARY FOR FUTURE USERS')

kar <- read.csv('Data/karyotype data from Boyes update direct copy 2011-06-04.csv')

kar$Genus <- as.character(kar$Genus)
kar$Genus[which(kar$Genus=='Rhysops')]<-'Argentinomyia'  ## Iron out nomenclature problems
 
kar$binomial<-paste(kar$Genus, kar$species, sep='_')  ### Create binomial name from genus/species columns
 
kar$binomial<-sub('_1$','',kar$binomial)
kar$binomial<-sub('_2$','',kar$binomial)
kar$binomial<-sub(' (nr)$','',kar$binomial)
kar$binomial[grep('Hiatomyia_plutonia',kar$binomial)][1]<-'Cheilosia_plutonia'
kar$binomial[grep('Cheilosia_ruffipes',kar$binomial)][1]<-'Cheilosia_rufipes'
 
## Create non-duplicated dataset
kar1 <- kar[!duplicated(kar$binomial),c('TCL','Npairs','tcl1','tcl2','ar2','tcl3','ar3','tcl4','ar4')]
row.names(kar1) <- kar$binomial[!duplicated(kar$binomial)]  

#kar1 <- subset(kar1, row.names(kar1)%in%phy$tip.label)
kar2<-apply(kar1,2,as.numeric)  ### Convert entire dataset to numeric values
row.names(kar2)<-row.names(kar1)

kar1 <- as.data.frame(na.omit(kar2))

