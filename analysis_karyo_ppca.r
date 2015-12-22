
### Run karyotype PPCA analysis

phy$node.label <- 1:phy$Nnode

##  kd is DATASET OF ALL SPP WITH ALL KARYOTYPE DATA ONLY - FOR PPCA PURPOSES - DO NOT USE FOR NPAIR ANALYSES
nrow(kd <- na.omit(subset(allkpdat, !is.na(TCL) & !is.na(Npair) & binomial %in% phy1$tip.label, select=c('binomial','TCL','Npair','tcl1','tcl2','tcl3','tcl4','ar2','ar3','ar4'))))#))))
# 303

Ntip(phyk <-drop.tip(compute.brlen(phy1, power=1), tip=which(!phy1$tip.label %in% kd$binomial))) ## compute branch lengths using Grafen's method
# 303
phyk$node.label <- 1:phyk$Nnode

ktreedat <- comparative.data(phyk, kd, binomial)

C<-vcv.phylo(ktreedat$phy)
C<-C[order(dimnames(C)[[1]]),order(dimnames(C)[[2]])]
kphy.dat <- as.matrix(ktreedat$data)
row.names(kphy.dat)<-row.names(ktreedat$data)
kpca.results<-phyl_pca(C, kphy.dat ,mode='corr')
ktreedat$data$Comp.1 <- as.vector(kpca.results$S[,1])
ktreedat$data$Comp.2 <- as.vector(kpca.results$S[,2])
ktreedat$data$Comp.3 <- as.vector(kpca.results$S[,3])
ktreedat$data$Comp.4 <- as.vector(kpca.results$S[,4])

### examine loadings for 1st 4 axes
kload <- kpca.results$L[,1:4]
row.names(kload)<-names(ktreedat$data[1:9])
kload
#     [,1]        [,2]        [,3]        [,4]
#     TCL   -0.32463248  0.50907728 -0.48888619 -0.26784116## PPCA 1 is almost perfectly correlated with Npair
#     Npair -0.93259802  0.04544417 -0.05797572 -0.09834670## Npair is closely negatively related to TCL2, 3 and 4 
#     tcl1   0.22631811 -0.68575969  0.24563895 -0.16690782
#     tcl2   0.83720912  0.18520421  0.07829307  0.01388209## PPCA2 is TCL and -tcl1
#     tcl3   0.85335240  0.20339893  0.08470288  0.08981079
#     tcl4   0.83471895  0.01068293 -0.21927353  0.18639664
#     ar2   -0.56349546 -0.21989584  0.32248979  0.24286087
#     ar3   -0.32322769  0.42147659  0.31516280  0.64474886
#     ar4   -0.08316971 -0.45237292 -0.72652232  0.47485087             

diag(kpca.results$Eval)/sum(diag(kpca.results$Eval))
#	[1] 0.39790565 0.13754611 0.12169518 0.09471527 0.09165409 0.07417733 0.04445459 0.02671898 0.01113281

