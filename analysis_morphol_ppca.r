
### Create morphological PPCA

#### First extract PPCA axes from morphological dataset 
phy1$node.label <- 1:phy1$Nnode
nrow(md <- subset(allkpdat, !is.na(WL) & binomial %in% phy1$tip.label, select=intersect(names(allkpdat),names(mph)[-ncol(mph)]))) 	
## select all columns of allkpdat that are in mph, except for the last one (binomial.old)
# 327
Ntip(phym <-drop.tip(compute.brlen(phy1, power=1), tip=which(!phy1$tip.label %in% md$binomial))) ## compute branch lengths using Grafen's method
# 327
phym$node.label <- 1:phym$Nnode  ### Tree must have node labels
mtreedat <- comparative.data(phym, md, binomial)  ## runs "na.omit" - make sure have removed "genus" and "species" columns which are all NA
C<-vcv.phylo(mtreedat$phy)  ## Run the PPCA
C<-C[order(dimnames(C)[[1]]),order(dimnames(C)[[2]])]
phy.dat <- as.matrix(mtreedat$data[,2:18])  ### Later, need to try other branch length transformations too
row.names(phy.dat)<-row.names(mtreedat$data)
mpca.results<-phyl_pca(C, phy.dat ,mode='corr')
mtreedat$data$Comp.1 <- as.vector(mpca.results$S[,1])
mtreedat$data$Comp.2 <- as.vector(mpca.results$S[,2])
mtreedat$data$Comp.3 <- as.vector(mpca.results$S[,3])
mtreedat$data$Comp.4 <- as.vector(mpca.results$S[,4])


#  Look at PPCA axes for morphology and karyotype

### examine loadings for 1st 4 axes
mload <- mpca.results$L[,1:4]
row.names(mload)<-names(mtreedat$data[2:18])
mload
# as expected PPCA 1 is size: large same-sign loadings for all vars
# PPCA2 is relative proboscis length: biggest loadings for +PS, -PL, -LABR, -PR, -FU
# PPCA3 is abdomen shape: biggest loadings for T2, T3 and T4
#       [,1]        [,2]        [,3]        [,4]
#       WL    -0.9379219  0.17088643 -0.12683145 -0.08184388
#       WW    -0.9363765  0.14373794 -0.04388413 -0.06374284
#       HW    -0.9640032  0.09224181 -0.02474899 -0.03054917
#       PL    -0.7885317 -0.54885130 -0.02477991  0.11268394
#       FU    -0.8230904 -0.48508996 -0.02069468  0.11292882
#       THW   -0.9682532  0.08163783 -0.05469598 -0.03789311
#       THL   -0.9578699  0.14955681 -0.14944462 -0.08348661
#       THH   -0.9687275  0.08285665 -0.11473618 -0.10602442
#       THVOL -0.9027466  0.05986960 -0.07989230 -0.16263719
#       T2    -0.9167298  0.05785613  0.36331735 -0.01085868
#       T3    -0.9155129  0.08792941  0.37979083 -0.02604234
#       T4    -0.8841434  0.11849227  0.41520942 -0.01187569
#       HTL   -0.8853842  0.15046677 -0.23065040 -0.21014265
#       LABR  -0.8185445 -0.55127196 -0.09101417  0.04590524
#       PR    -0.8313154 -0.50799822 -0.05503861  0.10810616
#       LL    -0.9307745  0.21090427 -0.02539568  0.20067141
#       PS    -0.6500317  0.61131108 -0.15369445  0.40984246

diag(mpca.results$Eval)/sum(diag(mpca.results$Eval))  ### Relative importance of each axis
#	 [1] 0.7933899371 0.0977827823 0.0354093044 0.0206102657 0.0169713490 0.0071686739 0.0066331194 0.0056697240 0.0036193575 0.0030079483 0.0026167689 0.0021735668
#	[13] 0.0017607919 0.0013077145 0.0009000419 0.0006164649 0.0003621897



