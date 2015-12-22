
### Import morphological dataset 

prob <- read.table('Data/females2.txt', sep='\t', header=T) ### Morphology dataset inc. proboscis characters

prob$binomial<-paste(prob$GENUS, prob$SPECIES, sep='_')

### prob contains some duplicates - for now, just choose the mean

## prob1 is non-phylo data (longer dataset but with nonphylo PCA)
prob1 <- aggregate(prob, list(prob$binomial), mean)[,c(1, 10:26)]
row.names(prob1) <- prob1$Group.1
prob1<-prob1[,-1]

### Not used - preparing for a non-phylo PCA analysis
# pcprob <- princomp(prob1)  ### NON-PHYLOGENETIC PCA
# prob.dat <- as.data.frame(pcprob$scores)
# prob.dat <- as.data.frame(prob2[,c('Comp.1','Comp.2','names')])
# prob.dat$TCL <- kar1$TCL[match(row.names(prob.dat), row.names(kar1))]
# prob.dat$Npair <- kar1$Npair[match(row.names(prob.dat), row.names(kar1))]
 

### prob2 is PHYLOGENETICALLY SYNCED AND WITH PHYLO PCA
phy.dat <- drop.tip(phy, tip=which(!phy$tip.label%in%row.names(prob1)))
prob2 <- as.matrix(subset(prob1, row.names(prob1)%in%phy.dat$tip.label))

### Create phylogenetic PCA for prob2
C<-vcv.phylo(phy.dat)
C<-C[order(dimnames(C)[[1]]),order(dimnames(C)[[2]])]
pca.results<-phyl_pca(C,prob2,mode='corr')
prob2 <- as.data.frame(prob2[order(row.names(prob2)),])
prob2$Comp.1 <- as.vector(pca.results$S[,1])  ### First PPCA axis
prob2$Comp.2 <- as.vector(pca.results$S[,2])  ### Second PPCA axis
prob2$names <- row.names(prob2)

##### Further manipulations from 2014.syrphid.analysis.r - not used as far as I can see
#   
#   prob.dat <- as.data.frame(pca.results$S[,1:4])  ## 2014: try using phylo PCA scores...
#   names(prob.dat)<-c('Comp.1','Comp.2','Comp.3','Comp.4')
#   
#   ### Add karyotype data to proboscis dataset
#   prob.dat$TCL <- kar1$TCL[match(row.names(prob.dat), row.names(kar1))]
#   prob.dat$Npair <- kar1$Npair[match(row.names(prob.dat), row.names(kar1))]
#   prob.dat$tcl1 <- kar1$tcl1[match(row.names(prob.dat), row.names(kar1))]
#   prob.dat$tcl2 <- kar1$tcl2[match(row.names(prob.dat), row.names(kar1))]
#   prob.dat$tcl3 <- kar1$tcl3[match(row.names(prob.dat), row.names(kar1))]
#   prob.dat$tcl4 <- kar1$tcl4[match(row.names(prob.dat), row.names(kar1))]
#   prob.dat$ar2 <- kar1$ar2[match(row.names(prob.dat), row.names(kar1))]
#   prob.dat$ar3 <- kar1$ar3[match(row.names(prob.dat), row.names(kar1))]
#   prob.dat$ar4 <- kar1$ar4[match(row.names(prob.dat), row.names(kar1))]
#   
#   dat.tree <- subset(prob.dat, !is.na(TCL) & row.names(prob.dat) %in% phy$tip.label) # & Comp.1 >= -100)  # remove the huge fly (Myathropa) 
#   tree.dat <- drop.tip(phy, tip=which(!phy$tip.label %in% row.names(dat.tree)))
#   tree.dat$node.label <- 1 : Nnode(tree.dat)  ### Tree needs unique node names
#   
#   dat.tree$names <- row.names(dat.tree) 
#   dat.tree$rcomp1 <- ((-(dat.tree$Comp.1 + abs(min(dat.tree$Comp.1)))->a) + abs(min(a))+1) 
#   dat.tree$rcomp2 <- ((-(dat.tree$Comp.2 + abs(min(dat.tree$Comp.2)))->a) + abs(min(a))+1) 
#   dat.tree$rcomp3 <- ((-(dat.tree$Comp.3 + abs(min(dat.tree$Comp.3)))->a) + abs(min(a))+1) 
#   dat.tree$rcomp4 <- ((-(dat.tree$Comp.4 + abs(min(dat.tree$Comp.4)))->a) + abs(min(a))+1) 
#   dat.tree$logrcomp1 <- log(dat.tree$rcomp1)
#   dat.tree$logrcomp2 <- log(dat.tree$rcomp2)
#   dat.tree$logrcomp3 <- log(dat.tree$rcomp3)
#   dat.tree$logrcomp4 <- log(dat.tree$rcomp4)
#   tree.dat.di <- multi2di(tree.dat) 
#   


