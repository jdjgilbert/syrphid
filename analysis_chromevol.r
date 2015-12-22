
##################################################################################################

Reconstruct chromosomal evolution



# Use ChromEvol to estimate most likely mode of chromosome evolution for syrphids based on observed phylogeny
# 
# Approach 1:
#   # Simulate observed genus distribution for N trees (growTree command in caper)
#   Over each one, simulate chromosome evolution with parameters obtained from ChromEvol
# For each of 1000 trees estimate the correspondence between genera and chromosome number
# 
# Approach 2:
#   Simulate most likely history of chromosome number evolution
# Count how often new genus evolved with each chromosome change
# Simulate histories with same parameters x 10,000
# Count again for each history
# Estimate probability of observed distribution of genera


Begin with approach 2

Create dataset of all spp with at least Npair data

knp IS DATASET OF ALL SPP IN TREE WITH AT LEAST NPAIR/TCL DATA -- USE IN CHROMOSOME RECON USING ChromEvol
ow(knp <- na.omit(subset(allkpdat, !is.na(TCL) & !is.na(Npair) & binomial %in% phy1$tip.label, select=c('binomial','TCL','Npair')))) 
488
ip(phyknp <-drop.tip(compute.brlen(phy1, power=1), tip=which(!phy1$tip.label %in% knp$binomial))) ## compute branch lengths using Grafen's method
488

Export chromosome numbers as FASTA format

<- paste(">", knp$binomial, sep='')
<- cbind(aa, knp$Npair)
rphidchr <- matrix(bb, ncol=2, byrow=F)
<- data.frame(a=matrix(t(syrphidchr), ncol=1))
ite.table(s1, file='~/Dropbox/2014/2014 Syrphidae/chromevol/syrphidchr.dat', sep='\t', quote=F, row.names=F, col.names=F)

Export tree as Newick
yknp$node.label <- NULL
ite.tree(phyknp, file='~/Dropbox/2014/2014 Syrphidae/chromevol/syrphidtree.phy')

Create chromevol_input.txt file with following content:
  ainType All_Models		## Tell program to test likelihood of all available models of chromosome evolution
ataFile syrphidchr.dat
reeFile syrphidtree.phy

Syntax in terminal (confirmed works)
pplications/chromEvol/chromEvol_source-2.0/chromEvol chromevol_input.txt 

Takes a while

# For this topology/branch lengths, lowest AIC is for LINEAR_RATE_NO_DUPL model (dAIC == 2 exactly (phew))

chr <- read.tree(   '../chromevol/RESULTS_ROOT_7_2014-09-18/LINEAR_RATE_NO_DUPL/mlAncestors.tree')
yeschr <- read.tree('../chromevol/RESULTS_ROOT_7_2014-09-18/LINEAR_RATE_NO_DUPL/posteriorAncestors.tree')

pdf('chromevol/RESULTS_ROOT_7/LINEAR_RATE_NO_DUPL/mlAncestors.pdf', height=50)  #
plot(compute.brlen(mlchr,1), cex=0.5, show.node.label=T)
dev.off()

pdf('chromevol/RESULTS_ROOT_7/LINEAR_RATE_NO_DUPL/posteriorAncestors.pdf', height=50) #
plot(compute.brlen(bayeschr,1), cex=0.5, show.node.label=T)
dev.off()


Extract character states at nodes and plot pretty phylogenies

Ancnodes <- as.numeric(sub('N.*-','',gsub('\\[|\\]','', mlchr$node.label)))

yesAncnodes <- sub('N[[:digit:]]*_','',gsub('\\[|\\]','', bayeschr$node.label))
ad(bayesNodes <- as.data.frame(do.call(rbind, lapply(lapply(sapply(bayesAncnodes, strsplit, "//"), strsplit, "-"), c, recursive=T))))
V1   V2 V3   V4
7-0.99//10-0    7 0.99 10    0
7-0.99//8-0     7 0.99  8    0
7-0.99//8-0.1   7 0.99  8    0
7-0.99//8-0.2   7 0.99  8    0
7-0.73//6-0.14  7 0.73  6 0.14
5-0.45//4-0.29  5 0.45  4 0.29

w.names(bayesNodes) <- paste('n', 1:nrow(bayesNodes), sep='')
mes(bayesNodes) <- c('State1','Prob1','State2','Prob2')

yesNodes$mlnode <- mlAncnodes  ## add ML nodes to dataset
ad(bayesNodes)
State1 Prob1 State2 Prob2 mlnode
n1      7  0.99     10     0      7
n2      7  0.99      8     0      7
n3      7  0.99      8     0      7
n4      7  0.99      8     0      7
n5      7  0.73      6  0.14      7
n6      5  0.45      4  0.29      5

yesNodes$Prob1 <- as.numeric(as.character(bayesNodes$Prob1))
yesNodes$Prob2 <- as.numeric(as.character(bayesNodes$Prob2))
yesNodes$prob3 <- with(bayesNodes, ((1 * (State1==3)) * Prob1) + ((1 * (State2==3)) * Prob2))
yesNodes$prob4 <- with(bayesNodes, ((1 * (State1==4)) * Prob1) + ((1 * (State2==4)) * Prob2))
yesNodes$prob5 <- with(bayesNodes, ((1 * (State1==5)) * Prob1) + ((1 * (State2==5)) * Prob2))
yesNodes$prob6 <- with(bayesNodes, ((1 * (State1==6)) * Prob1) + ((1 * (State2==6)) * Prob2))
yesNodes$prob7 <- with(bayesNodes, ((1 * (State1==7)) * Prob1) + ((1 * (State2==7)) * Prob2))
yesNodes$other <- 1 - apply(bayesNodes[,6:10],1,sum)


Plot ML phylogeny
pdf('../chromevol/RESULTS_ROOT_7_2014-09-18//LINEAR_RATE_NO_DUPL/Anctest.pdf', height=50)  
plot(compute.brlen(mlchr,1), cex=0.5, show.node.label=F, label.offset=1)
nodelabels(pch=21, cex=1, bg=c('green','red','orange','brown','darkgreen')[as.numeric(factor(bayesNodes$mlnode))])
tiplabels(pch=21, cex=1, bg=c('green','red','magenta','pink','orange','yellow','brown','darkgreen')[as.numeric(factor(knp$Npair[match(sub('-.*','',mlchr$tip.label),knp$binomial)]))])
dev.off()

Plot bayes phylogeny
pdf('../chromevol/RESULTS_ROOT_7_2014-09-18//LINEAR_RATE_NO_DUPL/Anctestbayes.pdf', height=50)  
plot(compute.brlen(mlchr,1), cex=0.5, show.node.label=F, label.offset=1)
nodelabels(cex=0.5, pie=as.matrix(bayesNodes[,6:11]), piecol=c('green','red','orange','brown','darkgreen','lightgrey'))
tiplabels(pch=21, cex=1, bg=c('green','red','magenta','pink','orange','yellow','brown','darkgreen')[as.numeric(factor(knp$Npair[match(sub('-.*','',mlchr$tip.label),knp$binomial)]))])
dev.off()

Plot ML phylogeny under DEMI assumption
pdf('../chromevol/RESULTS_ROOT_7_2014-09-18//LINEAR_RATE_DEMI/Anctest.pdf', height=50)  
plot(compute.brlen(mlchr,1), cex=0.5, show.node.label=F, label.offset=1)
nodelabels(pch=21, cex=1, bg=c('green','red','orange','brown','darkgreen')[as.numeric(factor(bayesNodes$mlnode))])
tiplabels(pch=21, cex=1, bg=c('green','red','magenta','pink','orange','yellow','brown','darkgreen')[as.numeric(factor(knp$Npair[match(sub('-.*','',mlchr$tip.label),knp$binomial)]))])
dev.off()

Plot bayes phylogeny under DEMI assumption
pdf('../chromevol/RESULTS_ROOT_7_2014-09-18//LINEAR_RATE_DEMI/Anctestbayes.pdf', height=50)  
plot(compute.brlen(mlchr,1), cex=0.5, show.node.label=F, label.offset=1)
nodelabels(cex=0.5, pie=as.matrix(bayesNodes[,6:11]), piecol=c('green','red','orange','brown','darkgreen','lightgrey'))
tiplabels(pch=21, cex=1, bg=c('green','red','magenta','pink','orange','yellow','brown','darkgreen')[as.numeric(factor(knp$Npair[match(sub('-.*','',mlchr$tip.label),knp$binomial)]))])
dev.off()


# create genus-level phylogeny for the figure

hy <- mlchr
hy$tip.label <- sub('_.*','',mphy$tip.label)
ip(mphy<-drop.tip(mphy, tip=which(mphy$tip.label[1:(Ntip(mphy)-1)] == mphy$tip.label[2:Ntip(mphy)] )))
111

hy <- bayeschr
hy$tip.label <- sub('_.*','',bphy$tip.label)
ip(bphy<-drop.tip(bphy, tip=which(bphy$tip.label[1:(Ntip(bphy)-1)] == bphy$tip.label[2:Ntip(bphy)] )))
111

ad(mphy$node.label)
] "[N1-7]"  "[N2-7]"  "[N5-7]"  "[N6-5]"  "[N7-4]"  "[N11-4]"
ad(bphy$node.label)
] "[N1_7-0.99//11-0]"   "[N2_7-0.99//8-0]"    "[N3_7-0.99//8-0]"    "[N4_7-0.99//8-0]"    "[N5_7-0.72//6-0.14]" "[N6_5-0.44//4-0.25]"

ngth(mlAncnodes <- as.numeric(sub('N.*-','',gsub('\\[|\\]','', mphy$node.label))))

ngth(bayesAncnodes <- sub('N[[:digit:]]*_','',gsub('\\[|\\]','', bphy$node.label)))
ad(bayesNodes <- as.data.frame(do.call(rbind, lapply(lapply(sapply(bayesAncnodes, strsplit, "//"), strsplit, "-"), c, recursive=T))))
V1   V2 V3   V4
7-0.99//11-0    7 0.99 11    0
7-0.99//8-0     7 0.99  8    0
7-0.99//8-0.1   7 0.99  8    0
7-0.99//8-0.2   7 0.99  8    0
7-0.72//6-0.14  7 0.72  6 0.14
5-0.44//4-0.25  5 0.44  4 0.25

w.names(bayesNodes) <- paste('n', 1:nrow(bayesNodes), sep='')
mes(bayesNodes) <- c('State1','Prob1','State2','Prob2')

yesNodes$mlnode <- 1 # mlAncnodes  ## add ML nodes to dataset
ad(bayesNodes)
State1 Prob1 State2 Prob2 mlnode
n1      7  0.99     10     0      7
n2      7  0.99      8     0      7
n3      7  0.99      8     0      7
n4      7  0.99      8     0      7
n5      7  0.73      6  0.14      7
n6      5  0.45      4  0.29      5

yesNodes$Prob1 <- as.numeric(as.character(bayesNodes$Prob1))
yesNodes$Prob2 <- as.numeric(as.character(bayesNodes$Prob2))
yesNodes$prob3 <- with(bayesNodes, ((1 * (State1==3)) * Prob1) + ((1 * (State2==3)) * Prob2))
yesNodes$prob4 <- with(bayesNodes, ((1 * (State1==4)) * Prob1) + ((1 * (State2==4)) * Prob2))
yesNodes$prob5 <- with(bayesNodes, ((1 * (State1==5)) * Prob1) + ((1 * (State2==5)) * Prob2))
yesNodes$prob6 <- with(bayesNodes, ((1 * (State1==6)) * Prob1) + ((1 * (State2==6)) * Prob2))
yesNodes$prob7 <- with(bayesNodes, ((1 * (State1==7)) * Prob1) + ((1 * (State2==7)) * Prob2))
yesNodes$other <- 1 - apply(bayesNodes[,6:10],1,sum)

script to get tip labels (from 2014.syrphid.graphics.1.r)
rnd$genus <- sub('_.*','',karnd$binomial)
ow(a <- with(karnd, ftable(genus,Npair)))
ad(a <- matrix(with(karnd, ftable(genus,sub('4.5','5', sub('5.5','6', sub('4.66666666666667','5', Npair))))), c(106,6)))
w.names(a)<-sort(unique(karnd$genus))
lnames(a)<-sort(unique(sub('4.5','5', sub('5.5','6', sub('4.66666666666667','5', karnd$Npair)))))
<- as.data.frame(a)
main <- apply(subset(a, select=c(1:5)), 1, function(x) names(a)[which.max(x)])
prop <- apply(subset(a, select=c(1:5)), 1, function(x) max(x)/sum(x))
occ <- as.vector(table(karnd$genus))
<- subset(a, select=c(1:5))
<- t(apply(b, 1, function(x) x/sum(x)))  ## matrix of proportions within genus represented by that chromosome count

hy <- read.tree(text = write.tree(bphy))  ## Weird fudge to get tip labels in the correct order in the R object!

brary(Hmisc)

pdf('../chromevol/RESULTS_ROOT_7_2014-09-18//LINEAR_RATE_NO_DUPL/Anctestbayesgenus.pdf', height=25)  
plot(compute.brlen(bphy,power=0.5), cex=1, show.node.label=F, label.offset=0.1)
nodelabels(cex=0.75, pie=as.matrix(bayesNodes[,6:11]), piecol=c('green','red','yellow','black','darkgreen','lightgrey'))
tiplabels(cex=0.75,pie=b[match(bphy$tip.label,row.names(b)),], piecol=c('green','red','yellow','black','darkgreen'))
#tiplabels(pch=21, cex=1, bg=c('green','red','magenta','pink','orange','yellow','brown','darkgreen')[as.numeric(factor(knp$Npair[match(sub('-.*','',mlchr$tip.label),knp$binomial)]))])

subplot(fun = {
  curve(9.08038+4.19813*x, from=3, to=7, lwd=1, lty=2, ylim=c(0, 40), yaxt='n', xaxt='n', ylab='Rate\n(arbitrary units)', xlab='Chromosomes',mgp=c(0.75, 0.5, 0), cex.lab=0.75, cex.axis=0.5);
  axis(1, at=c(3:7), tck=F, mgp=c(0,0,0), cex.axis=0.75);
  curve(17.7321-0.96324*x, from=3, to=7, add=T, lty=1, lwd=2);
  text(6,9.08038+4.19813*5+5, 'Chromosome loss', pos=2, cex=0.65);
  text(6,17.7321-0.96324*5+5, 'Chromosome gain', pos=2, cex=0.65);
}, x=0.05, y=105, size=c(1.5,1.5), hadj=0)

dev.off()



### Now have to ask the question: do new chromosomes == new genus?

## Create chromEvol simulation control file
_mainType mainSimulate		
_outDir RESULTS_SIMULATION/
  _treeFile syrphidtree.phy
_dataFile syrphidchr.txt
_freqFile rootfreq.txt
_lossConstR  1.52033   
_gainConstR  4.957043  
_lossLinearR 1.116383  
_gainLinearR -0.2711258

_simulationsIter 100
_maxChrNumForSimulations 17
_branchMul 1

## Parameters calculated from RESULTS_ROOT_NR/LINEAR_RATE_NO_DUPL/chromEvol.res
#	7.75654 *0.260245	#	LOSS_CONST      
#	27.3892 *0.260245	#	GAIN_CONST      
#	4.39315 *0.260245	#	LOSS_LINEAR     
#	-1.4911 *0.260245	#	GAIN_LINEAR     ## For some reason, says "linear gain paramater out of bounds"

## Parameters calculated from RESULTS_ROOT_7/LINEAR_RATE_NO_DUPL/chromEvol.res
5.84192 *0.260245	#	LOSS_CONST 
19.0476 *0.260245	#	GAIN_CONST      
4.28974 *0.260245	#	LOSS_LINEAR     
-1.04181 *0.260245	#	GAIN_LINEAR     
## gives the following parameters
_lossConstR  1.52033   
_gainConstR  4.957043  
_lossLinearR 1.116383  
_gainLinearR -0.2711258

## Parse output
## For some reason there is no terminal ";" on the trees output by the simulations

# Make sure trees are all identical
Nnode(t0 <- read.tree(text=paste(readLines('chromevol/RESULTS_SIMULATION/0/simTree.phr', 1),';', sep='')))
# 263
Nnode(t1 <- read.tree(text=paste(readLines('chromevol/RESULTS_SIMULATION/1/simTree.phr', 1),';', sep='')))
# 263
t0$tip.label <- 1:Ntip(t0); t1$tip.label <- 1:Ntip(t1)
t0$node.label <- 1:Nnode(t0); t1$node.label <- 1:Nnode(t1)
identical(t0,t1) ## TRUE

simtrees <- list()
for (i in 0:499) {
  tr <- try(read.tree(text=paste(readLines(paste('chromevol/RESULTS_SIMULATION/', i, '/simTree.phr', sep=''), 1),';', sep='')), silent=T)
  n <- data.frame(genus=sub('_.*','',tr$tip.label), sp=sub('.*?_','',sub('-.*','',tr$tip.label)), np=as.numeric(sub('.*-','',tr$tip.label)))
  ## Note: *? matches non-greedily!
  tr$tip.label <- paste(n$genus, n$sp, sep='_')
  simtrees[[i+1]] <- tr
  simtrees[[i+1]]$np <- n
  warning(i, immediate.=T)
}

res <- list()

for (i in 0:99) {
  tr <- read.tree(text=paste(readLines(paste('chromevol/RESULTS_SIMULATION/', i, '/simTree.phr', sep=''), 1),';', sep=''))
  n <- data.frame(genus=sub('_.*','',tr$tip.label), sp=sub('.*_','',sub('-.*','',tr$tip.label)), np=as.numeric(sub('.*-','',tr$tip.label)))
  head(a<-with(n, data.frame(table(genus,np))))
  head(a<-subset(a, Freq>0))
  head(b<-table(a$genus))
  res[[i+1]] <- as.vector(table(b))
}

max(unlist(lapply(res, length)))
res1 <- lapply(c(1:100), function(x) c(res[[x]], rep(0, 6-length(res[[x]]))))

res2 <- do.call(rbind, res1)
boxplot(res2)

## Compare with observed data
n <- data.frame(genus=sub('_.*','',knp$binomial), np=knp$Npair)
head(a<-with(n, data.frame(table(genus,np))))
head(a<-subset(a, Freq>0))
head(b<-table(a$genus))
table(b)

boxplot(t(data.frame(as.vector(table(b)))), add=T, border='red', col='red')
### Interesting -- appears to be a DEFICIT of genera with only 1 chromosome type....
### ... and more than expected genera with 2 types..
### Need to check assumptions etc. - branch lengths?


### Confirm this using outgroup - what is sister group to Syrphidae and what does their karyotype look like?
## tolweb.org --> Sister group is the Pipunculidae
## Google search for "Pipunculidae karyotype chromosome number" fruitless
### Pseudactaeon (Phoridae) has 2n = 6 (so n=3) - Chirino et al 2009, Genet Mol Biol
### Megaseila (Phoridae) has 2n = 6 (so n=3) - Wolf et al 1996, Hereditas

### Dad --> Microdon effectively an outgroup - some ppl consider a different family



## for context, other insects (from Wikipedia..)
## Honeybee n = 16 (8 for males)
## Myrmecia pilosula ants n=2 (1 for males)
## Drosophila n = 8
## All Culicidae n = 3 except Chagasia bathana n = 4
## Agriodiaetus butterfly 2n = 268 (so n = 134) one of highest among multicelled animals...



###### COMPARE RESULTS OF CHROMOSOME SIMULATION ON ALLOMETRY

simres <- list()

for (i in 1:500) {
  
  data <- simtrees[[i]]$np
  
  row.names(data) <- paste(data$genus, data$sp, sep='_')
  data$binomial <- paste(data$genus, data$sp, sep='_')
  data$Comp.2 <- mtreedat$data$Comp.2[match(row.names(data), row.names(mtreedat$data))]
  data$Comp.1 <- mtreedat$data$Comp.1[match(row.names(data), row.names(mtreedat$data))]
  data <- subset(data, !is.na(data$Comp.2))
  tree <- drop.tip(simtrees[[i]], tip=which(!simtrees[[i]]$tip.label%in%row.names(data)))
  cdata <- comparative.data(phy=tree,data=data,names='binomial')
  #m1 <- pgls(Comp.2 ~ Comp.1 * np, data=cdata, lambda='ML')
  #m0 <- pgls(Comp.2 ~ Comp.1 + np, data=cdata, lambda=m1$param[2])
  m1 <- gls(Comp.2 ~ Comp.1 * np, data=cdata$data, correlation=corPagel(1, cdata$phy),  method='ML')
  m0 <- gls(Comp.2 ~ Comp.1 + np, data=cdata$data, correlation=corPagel(1, cdata$phy),  method='ML')
  
  simres[[i]] <- m1
  simres[[i]]$m0 <- m0
  simres[[i]]$anova <- anova(m1, m0)
  simres[[i]]$lrt <- lrt(m1, m0)
  
  warning(i, immediate.=T)
  
}

simstats <- sapply(simres, function(x) x$anova$L.Ratio[2])
simstats <- sapply(simres, function(x) x$lrt$LR)


hist(simstats, breaks=20)
abline(v=quantile(simstats, 0.95))
abline(v=10.73, col='red')  ## LR test statistic from Npair model above


