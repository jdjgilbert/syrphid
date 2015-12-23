
  
## Produce big figure - genera tree with allometry plots next to it

## First, dataset with full data for Npair
#  use dataset "knp" as used in the chromevol analysis

nrow(knp <- na.omit(subset(allkpdat, !is.na(TCL) & !is.na(Npair) & binomial %in% phy1$tip.label, select=c('binomial','TCL','Npair')))) 
## 489 in 2014 annotations
## 488 when I do this in new Git project - trace this if necessary

Ntip(phyknp <-drop.tip(compute.brlen(phy1, power=1), tip=which(!phy1$tip.label %in% knp$binomial))) ## compute branch lengths using Grafen's method
## 489 in 2014 annotation
## 488 when I do this in new Git project - trace this if necessary
 

tree <- phyknp
tree$tip.label <- sub('_.*','',tree$tip.label)
Ntip(tree<-drop.tip(tree, tip=which(duplicated(tree$tip.label))))

## 106 (using phyknp)
# 81 (using mtreedat$phy)
# 64 (using kmtreedat$phy)
# 77 (using ktreedat$phy)


Ntip(tree1<-drop.tip(tree, tip=which(tree$tip.label[1:(Ntip(tree)-1)] == tree$tip.label[2:Ntip(tree)] )))
## 111
## 106 when I do this in new Git project - trace this if necessary

## get a from 2014.syrphid.analysis.X.r - table of chromosome counts per genus


### START here
pdf('Fig X - genera chromosome allometry 2014-12-01.pdf', width=12, height=25)

require(Hmisc)

karnd$genus <- sub('_.*','',karnd$binomial)
nrow(a <- with(karnd, ftable(genus,Npair)))
head(a <- matrix(with(karnd, ftable(genus,sub('4.5','5', sub('5.5','6', sub('4.66666666666667','5', Npair))))), c(106,6)))
row.names(a)<-sort(unique(karnd$genus))
colnames(a)<-sort(unique(sub('4.5','5', sub('5.5','6', sub('4.66666666666667','5', karnd$Npair)))))
a <- as.data.frame(a)

a$main <- apply(subset(a, select=c(1:5)), 1, function(x) names(a)[which.max(x)])
a$prop <- apply(subset(a, select=c(1:5)), 1, function(x) max(x)/sum(x))
a$occ <- as.vector(table(karnd$genus))

b <- subset(a, select=c(1:5))
b <- t(apply(b, 1, function(x) x/sum(x)))  ## matrix of proportions within genus represented by that chromosome count

#tree <- tree1
#tree <- read.tree(text = write.tree(tree))  ## Weird fudge to get tip labels in the correct order in the R object!
#plot(tree, label.offset=0.01)
#tiplabels(cex=0.25,pie=b[match(tree$tip.label,row.names(b)),], piecol=c('green','red','yellow','black','darkgreen'))

#	nf <- layout(matrix(c(1,2,1,2), 2, 2, byrow = TRUE), respect=T)
#	layout.show(nf)

genus3 <- row.names(subset(a, main==3 & prop>0.65))  ## Select spp with > two thirds of one chromosome type
genus4 <- row.names(subset(a, main==4 & prop>0.65))
genus5 <- row.names(subset(a, main==5 & prop>0.65))
genus6 <- row.names(subset(a, main==6 & prop>0.65))
genus7 <- row.names(subset(a, main==7 & prop>0.65))

genusn <- table(mtreedat$data$genus) ## Table of genera with morph data
genus4 <- genus4[which(genus4 %in% names(genusn[genusn>4]))]
genus5 <- genus5[which(genus5 %in% names(genusn[genusn>4]))]
genus6 <- genus6[which(genus6 %in% names(genusn[genusn>4]))]


par(mfrow=c(1,2),mar=c(0,0,0,0))
#plot(compute.brlen(tree, power=0.45), label.offset=0.02)
#tiplabels(cex=0.5,pie=b[match(tree$tip.label,row.names(b)),], piecol=c('green','black','skyblue','red','darkgreen'))
#aa <- sapply(tree$tip.label, strwidth, units='figure', font=2)
#points(1.1+aa, c(1:Ntip(tree)))

#m <- rbind(c(0, 0.5, 0, 1), c(0.5, 1, 0, 1))
#split.screen(m)
#screen(1)

CORRECT <- 1.12

######  NEW CODE FROM 2014.syrphid.analysis.10.r  (get bphy etc. from this code)
tree <- bphy
aa <- sapply(tree$tip.label, strwidth, units='figure', font=3)
plot(compute.brlen(tree,power=0.5), cex=1, show.node.label=F, label.offset=0.03)
nodelabels(cex=0.75, pie=as.matrix(bayesNodes[,6:11]), piecol=c('green','black','skyblue','red','darkgreen','lightgrey'))
tiplabels(cex=0.75,pie=b[match(tree$tip.label,row.names(b)),], piecol=c('green','black','skyblue','red','darkgreen'))
#tiplabels(pch=21, cex=1, bg=c('green','red','magenta','pink','orange','yellow','brown','darkgreen')[as.numeric(factor(knp$Npair[match(sub('-.*','',mlchr$tip.label),knp$binomial)]))])
#######


## Draw lines from each genus with data
yl <- c(genus4,genus5,genus6)
yy <- which(tree$tip.label %in% yl & !duplicated(tree$tip.label))
yc <- c(rep('black',length(genus4)),rep('skyblue',length(genus5)),rep('red',length(genus6)))
#lapply(c(1:length(yy)), function(x) lines(c(1.16+aa[yy[x]], 2+aa[yy[x]]), rep(yy[x], 2), col=yc[match(tree$tip.label[yy], yl)][x], lwd=2))
lapply(c(1:length(yy)), function(x) lines(c(aa[yy[x]]+CORRECT, 5+aa[yy[x]]), rep(yy[x], 2), col=yc[match(tree$tip.label[yy], yl)][x], lwd=2))
text(0.02,113, '(a)', pos=4, cex=1.75);

subplot(fun = {
  curve(9.08038+4.19813*x, from=3, to=7, lwd=1, lty=2, ylim=c(0, 40), yaxt='n', xaxt='n', ylab='Rate\n(arbitrary units)', xlab='Chromosomes',mgp=c(0.75, 0.5, 0), cex.lab=0.75, cex.axis=0.5);
  axis(1, at=c(3:7), tck=F, mgp=c(0,0,0), cex.axis=0.75);
  curve(17.7321-0.96324*x, from=3, to=7, add=T, lty=1, lwd=2);
  text(6,9.08038+4.19813*5+5, 'Chromosome loss', pos=2, cex=0.65, font=2);
  text(6,17.7321-0.96324*5+5, 'Chromosome gain', pos=2, cex=0.65, font=2);
  text(2.75,5, '(b)', pos=4, cex=1.75);
}, x=0.07, y=105, size=c(1.5,1.5), hadj=0)


#screen(2)

## move to the next plot window but plot the tree in it, so the XY coords of both windows match 
par(mfg=c(1, 2))
plot(tree, plot=F)
#plot.new()
#lapply(c(1:Ntip(tree)), function(x) lines(c(-0.2, aa[x]), rep(x, 2)))

LEFT <- -0.0135
L4 <- 0.05
L5 <- L4 + 0.07
R5 <- L5 + 0.06
L6 <- R5 + 0.07
R6 <- L6 + 0.06

SUBPLOT_SIZE <- 0.75

## Draw lines to the relevant genera on the tree
## 4 chromosomes
##	yPipiza <- grep("Archimicrodon", tree$tip.label)
yNeocnemodon <- grep("Melanostoma", tree$tip.label)
ySphaerophoria <- grep("Ocyptamus", tree$tip.label)[1]
yPlatycheirus <- grep("Allograpta", tree$tip.label)[2]
ySyrphus <- grep("Austrosyrphus", tree$tip.label)
yParasyrphus <- grep("Betasyrphus", tree$tip.label)
yEupeodes <- grep("Ischiodon", tree$tip.label)
yTemnostoma <- grep("Criorhina", tree$tip.label)

##	lines(c(L4, LEFT), c(yPipiza, grep('Pipiza', tree$tip.label)), lwd=2)
lines(c(L4, LEFT), c(yNeocnemodon, grep('Neocnemodon', tree$tip.label)), lwd=2)
lines(c(L4, LEFT), c(yPlatycheirus, grep('Platycheirus', tree$tip.label)), lwd=2)
lines(c(L4, LEFT), c(ySphaerophoria, grep('Sphaerophoria', tree$tip.label)), lwd=2)
lines(c(L4, LEFT), c(yParasyrphus, grep('Parasyrphus', tree$tip.label)), lwd=2)
lines(c(L4, LEFT), c(ySyrphus, grep('Syrphus', tree$tip.label)), lwd=2)
lines(c(L4, LEFT), c(yEupeodes, grep('Eupeodes', tree$tip.label)), lwd=2)
lines(c(L4, LEFT), c(yTemnostoma, grep('Temnostoma', tree$tip.label)), lwd=2)
#Eupeodes
#Neocnemodon
#Parasyrphus

## 5 Chromosomes
##	yToxomerus <- grep('Ocyptamus', tree$tip.label)[3]
yXylota <- grep('Tropidia', tree$tip.label)
##	yEpistrophe <- grep('Xanthogramma', tree$tip.label)
##	yIschyrosyrphus <- grep('Pseudodoros', tree$tip.label)
##	##	yMelangyna <- grep('Ischyrosyrphus', tree$tip.label)
##	ySpilomyia <- grep('Milesia', tree$tip.label)
yBrachyopa <- grep('Volucella', tree$tip.label)
yChalcosyrphus <- grep('Hadromyia', tree$tip.label)
yChrysotoxum <- grep('Eupeodes', tree$tip.label)
yDasysyrphus <- grep('Didea', tree$tip.label)
##	lines(c(L5, LEFT), c(yToxomerus, grep('Toxomerus', tree$tip.label)), col='skyblue', lwd=2)
lines(c(L5, LEFT), c(yXylota, grep('Xylota', tree$tip.label)), col='skyblue', lwd=2)
##	lines(c(R5, LEFT), c(yEpistrophe, grep('Epistrophe$', tree$tip.label)), col='skyblue', lwd=2)
##	lines(c(L5, LEFT), c(yIschyrosyrphus, grep('Ischyrosyrphus', tree$tip.label)), col='skyblue', lwd=2)
##	lines(c(L5, LEFT), c(yMelangyna, grep('Melangyna', tree$tip.label)), col='skyblue', lwd=2)
##	lines(c(L5, LEFT), c(ySpilomyia, grep('Spilomyia', tree$tip.label)), col='skyblue', lwd=2)
lines(c(L5, LEFT), c(yBrachyopa, grep('Brachyopa', tree$tip.label)), col='skyblue', lwd=2)
lines(c(R5, LEFT), c(yChalcosyrphus, grep('Chalcosyrphus', tree$tip.label)), col='skyblue', lwd=2)
lines(c(R5, LEFT), c(yChrysotoxum, grep('Chrysotoxum', tree$tip.label)), col='skyblue', lwd=2)
lines(c(L5, LEFT), c(yDasysyrphus, grep('Dasysyrphus', tree$tip.label)), col='skyblue', lwd=2)
##	#Toxomerus
#Xylota
##	#Epistrophe
##	#Ischyrosyrphus
##	#Melangyna
##	#Spilomyia
#Brachyopa
#Chalcosyrphus
#Chrysotoxum
#Dasysyrphus

## 6 Chromosomes
##	yCheilosia <- grep('Merodon', tree$tip.label)
yHelophilus <- grep('Eurimyia', tree$tip.label)
ySphegina <- grep('Neoascia', tree$tip.label)
yBlera <- grep('Eristalinus', tree$tip.label)
yChrysogaster <- grep('Brachyopa', tree$tip.label)
yEoseristalis <- grep('Pilinascia', tree$tip.label)
##	yCopestylum <- grep('Ferdinandea', tree$tip.label)
yCriorhina <- grep('Somula', tree$tip.label)
ySericomyia <- grep('Simoides', tree$tip.label)
ySphecomyia <- grep('Myolepta', tree$tip.label)
yMallota <- grep('Myathropa', tree$tip.label)

##	lines(c(L6, LEFT), c(yCheilosia, grep('Cheilosia', tree$tip.label)[1]), col='red', lwd=2)
lines(c(L6, LEFT), c(yHelophilus, grep('Helophilus', tree$tip.label)), col='red', lwd=2)
lines(c(L6, LEFT), c(ySphegina, grep('Sphegina', tree$tip.label)), col='red', lwd=2)
lines(c(L6, LEFT), c(yBlera, grep('Blera', tree$tip.label)), col='red', lwd=2)
lines(c(L6, LEFT), c(yMallota, grep('Mallota', tree$tip.label)), col='red', lwd=2)
#lines(c(L6, LEFT), c(yChrysogaster, grep('Chrysogaster', tree$tip.label)), col='red', lwd=2)
lines(c(R6, LEFT), c(yEoseristalis, grep('Eoseristalis', tree$tip.label)), col='red', lwd=2)
##	lines(c(L6, LEFT), c(yCopestylum, grep('Copestylum', tree$tip.label)), col='red', lwd=2)
lines(c(L6, LEFT), c(yCriorhina, grep('Criorhina', tree$tip.label)), col='red', lwd=2)
lines(c(R6, LEFT), c(ySericomyia, grep('Sericomyia', tree$tip.label)), col='red', lwd=2)
lines(c(R6, LEFT), c(ySphecomyia, grep('Sphecomyia', tree$tip.label)), col='red', lwd=2)
#Helophilus
#Sphegina
#Blera
#Chrysogaster
#Eosristalis

## Produce subplots

par(mfg=c(1, 2))
## draw legend and schematic plot area
legend(0.20, 20, fill=c('green','black','skyblue','red','darkgreen'), legend=c('3 chromosomes','4 chromosomes','5 chromosomes','6 chromosomes','7 chromosomes'), bty='n', inset=0.15, title='Key')
subplot(fun = {
  plot(c(1:10)~c(1:10), type='n', tck=0.05, xaxt='n', yaxt='n', ylab='Shape', xlab='Size', mgp=c(0, 0.25, 0));
  axis(1, at=c(1:10), tck=0.05, labels=F);
  axis(2, at=c(1:10), tck=0.05, labels=F)
}, x=0.225, y=10, size=c(1,1), hadj=0)


par(mfg=c(1, 2))
##	subplot(fun = {
##	with(mtreedat$data[which(mtreedat$data$genus %in% 'Pipiza'),], plot(Comp.2 ~ Comp.1, pch=16, cex=0.5, tck=F,xlab='', ylab='', xaxt='n', yaxt='n', ylim=range(mtreedat$data$Comp.2[mtreedat$data$genus%in%genus4]), xlim=range(mtreedat$data$Comp.1[mtreedat$data$genus%in%genus4])));
##	rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white");
##	axis(1, mgp=c(0, -0.25, 0), labels=T, cex.axis=0.5, at=c(-0.4, 0.2), tck=0.05);
##	axis(2, mgp=c(0, 0.25, 0), cex.axis=0.5, at=c(-0.1, 0.2), tck=0.05, las=1);
##	with(mtreedat$data[which(mtreedat$data$genus %in% 'Pipiza'),], points(Comp.2 ~ Comp.1, pch=16, cex=0.5))
##	}, x=L4, y=yPipiza, size=c(SUBPLOT_SIZE,SUBPLOT_SIZE), hadj=0)

subplot(fun = {
  with(mtreedat$data[which(mtreedat$data$genus %in% 'Neocnemodon'),], plot(Comp.2 ~ Comp.1, pch=16, cex=0.5, tck=F,xlab='', ylab='', xaxt='n', yaxt='n', ylim=range(mtreedat$data$Comp.2[mtreedat$data$genus%in%genus4]), xlim=range(mtreedat$data$Comp.1[mtreedat$data$genus%in%genus4])));
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white");
  with(mtreedat$data[which(mtreedat$data$genus %in% 'Neocnemodon'),], points(Comp.2 ~ Comp.1, pch=16, cex=0.5))
}, x=L4, y=yNeocnemodon, size=c(SUBPLOT_SIZE,SUBPLOT_SIZE), hadj=0)

subplot(fun = {
  with(mtreedat$data[which(mtreedat$data$genus %in% 'Platycheirus'),], plot(Comp.2 ~ Comp.1, pch=16, cex=0.5, tck=F,xlab='', ylab='', xaxt='n', yaxt='n', ylim=range(mtreedat$data$Comp.2[mtreedat$data$genus%in%genus4]), xlim=range(mtreedat$data$Comp.1[mtreedat$data$genus%in%genus4])));
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white");
  with(mtreedat$data[which(mtreedat$data$genus %in% 'Platycheirus'),], points(Comp.2 ~ Comp.1, pch=16, cex=0.5))
}, x=L4, y=yPlatycheirus, size=c(SUBPLOT_SIZE,SUBPLOT_SIZE), hadj=0)

subplot(fun = {
  with(mtreedat$data[which(mtreedat$data$genus %in% 'Sphaerophoria'),], plot(Comp.2 ~ Comp.1, pch=16, cex=0.5, tck=F,xlab='', ylab='', xaxt='n', yaxt='n', ylim=range(mtreedat$data$Comp.2[mtreedat$data$genus%in%genus4]), xlim=range(mtreedat$data$Comp.1[mtreedat$data$genus%in%genus4])));
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white");
  with(mtreedat$data[which(mtreedat$data$genus %in% 'Sphaerophoria'),], points(Comp.2 ~ Comp.1, pch=16, cex=0.5))
}, x=L4, y=ySphaerophoria, size=c(SUBPLOT_SIZE,SUBPLOT_SIZE), hadj=0)

subplot(fun = {
  with(mtreedat$data[which(mtreedat$data$genus %in% 'Parasyrphus'),], plot(Comp.2 ~ Comp.1, pch=16, cex=0.5, tck=F,xlab='', ylab='', xaxt='n', yaxt='n', ylim=range(mtreedat$data$Comp.2[mtreedat$data$genus%in%genus4]), xlim=range(mtreedat$data$Comp.1[mtreedat$data$genus%in%genus4])));
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white");
  with(mtreedat$data[which(mtreedat$data$genus %in% 'Parasyrphus'),], points(Comp.2 ~ Comp.1, pch=16, cex=0.5))
}, x=L4, y=yParasyrphus, size=c(SUBPLOT_SIZE,SUBPLOT_SIZE), hadj=0)

subplot(fun = {
  with(mtreedat$data[which(mtreedat$data$genus %in% 'Syrphus'),], plot(Comp.2 ~ Comp.1, pch=16, cex=0.5, tck=F,xlab='', ylab='', xaxt='n', yaxt='n', ylim=range(mtreedat$data$Comp.2[mtreedat$data$genus%in%genus4]), xlim=range(mtreedat$data$Comp.1[mtreedat$data$genus%in%genus4])));
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white");
  with(mtreedat$data[which(mtreedat$data$genus %in% 'Syrphus'),], points(Comp.2 ~ Comp.1, pch=16, cex=0.5))
}, x=L4, y=ySyrphus, size=c(SUBPLOT_SIZE,SUBPLOT_SIZE), hadj=0)

subplot(fun = {
  with(mtreedat$data[which(mtreedat$data$genus %in% 'Eupeodes'),], plot(Comp.2 ~ Comp.1, pch=16, cex=0.5, tck=F,xlab='', ylab='', xaxt='n', yaxt='n', ylim=range(mtreedat$data$Comp.2[mtreedat$data$genus%in%genus4]), xlim=range(mtreedat$data$Comp.1[mtreedat$data$genus%in%genus4])));
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white");
  with(mtreedat$data[which(mtreedat$data$genus %in% 'Eupeodes'),], points(Comp.2 ~ Comp.1, pch=16, cex=0.5))
}, x=L4, y=yEupeodes, size=c(SUBPLOT_SIZE,SUBPLOT_SIZE), hadj=0)

subplot(fun = {
  with(mtreedat$data[which(mtreedat$data$genus %in% 'Temnostoma'),], plot(Comp.2 ~ Comp.1, pch=16, cex=0.5, tck=F,xlab='', ylab='', xaxt='n', yaxt='n', ylim=range(mtreedat$data$Comp.2[mtreedat$data$genus%in%genus4]), xlim=range(mtreedat$data$Comp.1[mtreedat$data$genus%in%genus4])));
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white");
  with(mtreedat$data[which(mtreedat$data$genus %in% 'Temnostoma'),], points(Comp.2 ~ Comp.1, pch=16, cex=0.5))
}, x=L4, y=yTemnostoma, size=c(SUBPLOT_SIZE,SUBPLOT_SIZE), hadj=0)

## 5 Chromosomes

##	subplot(fun = {
##	with(mtreedat$data[which(mtreedat$data$genus %in% 'Toxomerus'),], plot(Comp.2 ~ Comp.1, pch=16, cex=0.5, col = 'skyblue', tck=F,xlab='', ylab='', xaxt='n', yaxt='n', ylim=range(mtreedat$data$Comp.2[mtreedat$data$genus%in%genus5]), xlim=range(mtreedat$data$Comp.1[mtreedat$data$genus%in%genus5])));
##	rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white");
##	axis(1, mgp=c(0, -0.25, 0), labels=T, cex.axis=0.5, at=c(-0.5, 0.3), tck=0.05);
##	axis(2, mgp=c(0, 0.25, 0), cex.axis=0.5, at=c(-0.1, 0.1), tck=0.05, las=1);
##	with(mtreedat$data[which(mtreedat$data$genus %in% 'Toxomerus'),], points(Comp.2 ~ Comp.1, pch=16, cex=0.5, col = 'skyblue'))
##	}, x=L5, y=yToxomerus, size=c(SUBPLOT_SIZE,SUBPLOT_SIZE), hadj=0)

subplot(fun = {
  with(mtreedat$data[which(mtreedat$data$genus %in% 'Xylota'),], plot(Comp.2 ~ Comp.1, pch=16, cex=0.5, col = 'skyblue', tck=F,xlab='', ylab='', xaxt='n', yaxt='n', ylim=range(mtreedat$data$Comp.2[mtreedat$data$genus%in%genus5]), xlim=range(mtreedat$data$Comp.1[mtreedat$data$genus%in%genus5])));
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white");
  with(mtreedat$data[which(mtreedat$data$genus %in% 'Xylota'),], points(Comp.2 ~ Comp.1, pch=16, cex=0.5, col = 'skyblue'))
}, x=L5, y=yXylota, size=c(SUBPLOT_SIZE,SUBPLOT_SIZE), hadj=0)


##	subplot(fun = {
##	with(mtreedat$data[which(mtreedat$data$genus %in% 'Epistrophe'),], plot(Comp.2 ~ Comp.1, pch=16, cex=0.5, col = 'skyblue', tck=F,xlab='', ylab='', xaxt='n', yaxt='n', ylim=range(mtreedat$data$Comp.2[mtreedat$data$genus%in%genus5]), xlim=range(mtreedat$data$Comp.1[mtreedat$data$genus%in%genus5])));
##	rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white");
##	with(mtreedat$data[which(mtreedat$data$genus %in% 'Epistrophe'),], points(Comp.2 ~ Comp.1, pch=16, cex=0.5, col = 'skyblue'))
##	}, x=R5, y=yEpistrophe, size=c(SUBPLOT_SIZE,SUBPLOT_SIZE), hadj=0)


##	subplot(fun = {
##	with(mtreedat$data[which(mtreedat$data$genus %in% 'Ischyrosyrphus'),], plot(Comp.2 ~ Comp.1, pch=16, cex=0.5, col = 'skyblue', tck=F,xlab='', ylab='', xaxt='n', yaxt='n', ylim=range(mtreedat$data$Comp.2[mtreedat$data$genus%in%genus5]), xlim=range(mtreedat$data$Comp.1[mtreedat$data$genus%in%genus5])));
##	rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white");
##	with(mtreedat$data[which(mtreedat$data$genus %in% 'Ischyrosyrphus'),], points(Comp.2 ~ Comp.1, pch=16, cex=0.5, col = 'skyblue'))
##	}, x=L5, y=yIschyrosyrphus, size=c(SUBPLOT_SIZE,SUBPLOT_SIZE), hadj=0)


##	subplot(fun = {
##	with(mtreedat$data[which(mtreedat$data$genus %in% 'Melangyna'),], plot(Comp.2 ~ Comp.1, pch=16, cex=0.5, col = 'skyblue', tck=F,xlab='', ylab='', xaxt='n', yaxt='n', ylim=range(mtreedat$data$Comp.2[mtreedat$data$genus%in%genus5]), xlim=range(mtreedat$data$Comp.1[mtreedat$data$genus%in%genus5])));
##	rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white");
##	with(mtreedat$data[which(mtreedat$data$genus %in% 'Melangyna'),], points(Comp.2 ~ Comp.1, pch=16, cex=0.5, col = 'skyblue'))
##	}, x=L5, y=yMelangyna, size=c(SUBPLOT_SIZE,SUBPLOT_SIZE), hadj=0)


#subplot(fun = {
##	with(mtreedat$data[which(mtreedat$data$genus %in% 'Spilomyia'),], plot(Comp.2 ~ Comp.1, pch=16, cex=0.5, col = 'skyblue', tck=F,xlab='', ylab='', xaxt='n', yaxt='n', ylim=range(mtreedat$data$Comp.2[mtreedat$data$genus%in%genus5]), xlim=range(mtreedat$data$Comp.1[mtreedat$data$genus%in%genus5])));
#	rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white");
##	with(mtreedat$data[which(mtreedat$data$genus %in% 'Spilomyia'),], points(Comp.2 ~ Comp.1, pch=16, cex=0.5, col = 'skyblue'))
##	}, x=L5, y=ySpilomyia, size=c(SUBPLOT_SIZE,SUBPLOT_SIZE), hadj=0)


subplot(fun = {
  with(mtreedat$data[which(mtreedat$data$genus %in% 'Brachyopa'),], plot(Comp.2 ~ Comp.1, pch=16, cex=0.5, col = 'skyblue', tck=F,xlab='', ylab='', xaxt='n', yaxt='n', ylim=range(mtreedat$data$Comp.2[mtreedat$data$genus%in%genus5]), xlim=range(mtreedat$data$Comp.1[mtreedat$data$genus%in%genus5])));
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white");
  with(mtreedat$data[which(mtreedat$data$genus %in% 'Brachyopa'),], points(Comp.2 ~ Comp.1, pch=16, cex=0.5, col = 'skyblue'))
}, x=L5, y=yBrachyopa, size=c(SUBPLOT_SIZE,SUBPLOT_SIZE), hadj=0)


subplot(fun = {
  with(mtreedat$data[which(mtreedat$data$genus %in% 'Chalcosyrphus'),], plot(Comp.2 ~ Comp.1, pch=16, cex=0.5, col = 'skyblue', tck=F,xlab='', ylab='', xaxt='n', yaxt='n', ylim=range(mtreedat$data$Comp.2[mtreedat$data$genus%in%genus5]), xlim=range(mtreedat$data$Comp.1[mtreedat$data$genus%in%genus5])));
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white");
  with(mtreedat$data[which(mtreedat$data$genus %in% 'Chalcosyrphus'),], points(Comp.2 ~ Comp.1, pch=16, cex=0.5, col = 'skyblue'))
}, x=R5, y=yChalcosyrphus, size=c(SUBPLOT_SIZE,SUBPLOT_SIZE), hadj=0)


subplot(fun = {
  with(mtreedat$data[which(mtreedat$data$genus %in% 'Chrysotoxum'),], plot(Comp.2 ~ Comp.1, pch=16, cex=0.5, col = 'skyblue', tck=F,xlab='', ylab='', xaxt='n', yaxt='n', ylim=range(mtreedat$data$Comp.2[mtreedat$data$genus%in%genus5]), xlim=range(mtreedat$data$Comp.1[mtreedat$data$genus%in%genus5])));
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white");
  with(mtreedat$data[which(mtreedat$data$genus %in% 'Chrysotoxum'),], points(Comp.2 ~ Comp.1, pch=16, cex=0.5, col = 'skyblue'))
}, x=R5, y=yChrysotoxum, size=c(SUBPLOT_SIZE,SUBPLOT_SIZE), hadj=0)


subplot(fun = {
  with(mtreedat$data[which(mtreedat$data$genus %in% 'Dasysyrphus'),], plot(Comp.2 ~ Comp.1, pch=16, cex=0.5, col = 'skyblue', tck=F,xlab='', ylab='', xaxt='n', yaxt='n', ylim=range(mtreedat$data$Comp.2[mtreedat$data$genus%in%genus5]), xlim=range(mtreedat$data$Comp.1[mtreedat$data$genus%in%genus5])));
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white");
  with(mtreedat$data[which(mtreedat$data$genus %in% 'Dasysyrphus'),], points(Comp.2 ~ Comp.1, pch=16, cex=0.5, col = 'skyblue'))
}, x=L5, y=yDasysyrphus, size=c(SUBPLOT_SIZE,SUBPLOT_SIZE), hadj=0)


## 6 Chromosomes

##	subplot(fun = {
##	with(mtreedat$data[which(mtreedat$data$genus %in% 'Cheilosia'),], plot(Comp.2 ~ Comp.1, pch=16, cex=0.5, col='red', tck=F,xlab='', ylab='', xaxt='n', yaxt='n', ylim=range(mtreedat$data$Comp.2[mtreedat$data$genus%in%genus6]), xlim=range(mtreedat$data$Comp.1[mtreedat$data$genus%in%genus6])));
##	axis(1, mgp=c(0, -0.25, 0), labels=T, cex.axis=0.5, at=c(-1, 0.2), tck=0.05);
##	axis(2, mgp=c(0, 0.25, 0), cex.axis=0.5, at=c(-0.2, 0), tck=0.05, las=1);
##	rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white");
##	with(mtreedat$data[which(mtreedat$data$genus %in% 'Cheilosia'),], points(Comp.2 ~ Comp.1, pch=16, cex=0.5, col='red'))
##	}, x=L6, y=yCheilosia, size=c(SUBPLOT_SIZE,SUBPLOT_SIZE), hadj=0)


subplot(fun = {
  with(mtreedat$data[which(mtreedat$data$genus %in% 'Helophilus'),], plot(Comp.2 ~ Comp.1, pch=16, cex=0.5, col='red', tck=F,xlab='', ylab='', xaxt='n', yaxt='n', ylim=range(mtreedat$data$Comp.2[mtreedat$data$genus%in%genus6]), xlim=range(mtreedat$data$Comp.1[mtreedat$data$genus%in%genus6])));
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white");
  with(mtreedat$data[which(mtreedat$data$genus %in% 'Helophilus'),], points(Comp.2 ~ Comp.1, pch=16, cex=0.5, col='red'))
}, x=L6, y=yHelophilus, size=c(SUBPLOT_SIZE,SUBPLOT_SIZE), hadj=0)

subplot(fun = {
  with(mtreedat$data[which(mtreedat$data$genus %in% 'Sphegina'),], plot(Comp.2 ~ Comp.1, pch=16, cex=0.5, col='red', tck=F,xlab='', ylab='', xaxt='n', yaxt='n', ylim=range(mtreedat$data$Comp.2[mtreedat$data$genus%in%genus6]), xlim=range(mtreedat$data$Comp.1[mtreedat$data$genus%in%genus6])));
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white");
  with(mtreedat$data[which(mtreedat$data$genus %in% 'Sphegina'),], points(Comp.2 ~ Comp.1, pch=16, cex=0.5, col='red'))
}, x=L6, y=ySphegina, size=c(SUBPLOT_SIZE,SUBPLOT_SIZE), hadj=0)

subplot(fun = {
  with(mtreedat$data[which(mtreedat$data$genus %in% 'Blera'),], plot(Comp.2 ~ Comp.1, pch=16, cex=0.5, col='red', tck=F,xlab='', ylab='', xaxt='n', yaxt='n', ylim=range(mtreedat$data$Comp.2[mtreedat$data$genus%in%genus6]), xlim=range(mtreedat$data$Comp.1[mtreedat$data$genus%in%genus6])));
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white");
  with(mtreedat$data[which(mtreedat$data$genus %in% 'Blera'),], points(Comp.2 ~ Comp.1, pch=16, cex=0.5, col='red'))
}, x=L6, y=yBlera, size=c(SUBPLOT_SIZE,SUBPLOT_SIZE), hadj=0)

#subplot(fun = {
#	with(mtreedat$data[which(mtreedat$data$genus %in% 'Chrysogaster'),], plot(Comp.2 ~ Comp.1, pch=16, cex=0.5, col='red', tck=F,xlab='', ylab='', xaxt='n', yaxt='n', ylim=range(mtreedat$data$Comp.2[mtreedat$data$genus%in%genus6]), xlim=range(mtreedat$data$Comp.1[mtreedat$data$genus%in%genus6])));
#	rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white");
#	with(mtreedat$data[which(mtreedat$data$genus %in% 'Chrysogaster'),], points(Comp.2 ~ Comp.1, pch=16, cex=0.5, col='red'))
#	}, x=L6, y=yChrysogaster, size=c(SUBPLOT_SIZE,SUBPLOT_SIZE), hadj=0)

subplot(fun = {
  with(mtreedat$data[which(mtreedat$data$genus %in% 'Eoseristalis'),], plot(Comp.2 ~ Comp.1, pch=16, cex=0.5, col='red', tck=F,xlab='', ylab='', xaxt='n', yaxt='n', ylim=range(mtreedat$data$Comp.2[mtreedat$data$genus%in%genus6]), xlim=range(mtreedat$data$Comp.1[mtreedat$data$genus%in%genus6])));
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white");
  with(mtreedat$data[which(mtreedat$data$genus %in% 'Eoseristalis'),], points(Comp.2 ~ Comp.1, pch=16, cex=0.5, col='red'))
}, x=R6, y=yEoseristalis, size=c(SUBPLOT_SIZE,SUBPLOT_SIZE), hadj=0)

##subplot(fun = {
##	with(mtreedat$data[which(mtreedat$data$genus %in% 'Copestylum'),], plot(Comp.2 ~ Comp.1, pch=16, cex=0.5, col='red', tck=F,xlab='', ylab='', xaxt='n', yaxt='n', ylim=range(mtreedat$data$Comp.2[mtreedat$data$genus%in%genus6]), xlim=range(mtreedat$data$Comp.1[mtreedat$data$genus%in%genus6])));
##rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white");
##	with(mtreedat$data[which(mtreedat$data$genus %in% 'Copestylum'),], points(Comp.2 ~ Comp.1, pch=16, cex=0.5, col='red'))
##	}, x=L6, y=yCopestylum, size=c(SUBPLOT_SIZE,SUBPLOT_SIZE), hadj=0)

subplot(fun = {
  with(mtreedat$data[which(mtreedat$data$genus %in% 'Criorhina'),], plot(Comp.2 ~ Comp.1, pch=16, cex=0.5, col='red', tck=F,xlab='', ylab='', xaxt='n', yaxt='n', ylim=range(mtreedat$data$Comp.2[mtreedat$data$genus%in%genus6]), xlim=range(mtreedat$data$Comp.1[mtreedat$data$genus%in%genus6])));
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white");
  with(mtreedat$data[which(mtreedat$data$genus %in% 'Criorhina'),], points(Comp.2 ~ Comp.1, pch=16, cex=0.5, col='red'))
}, x=L6, y=yCriorhina, size=c(SUBPLOT_SIZE,SUBPLOT_SIZE), hadj=0)

subplot(fun = {
  with(mtreedat$data[which(mtreedat$data$genus %in% 'Mallota'),], plot(Comp.2 ~ Comp.1, pch=16, cex=0.5, col='red', tck=F,xlab='', ylab='', xaxt='n', yaxt='n', ylim=range(mtreedat$data$Comp.2[mtreedat$data$genus%in%genus6]), xlim=range(mtreedat$data$Comp.1[mtreedat$data$genus%in%genus6])));
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white");
  with(mtreedat$data[which(mtreedat$data$genus %in% 'Mallota'),], points(Comp.2 ~ Comp.1, pch=16, cex=0.5, col='red'))
}, x=L6, y=yMallota, size=c(SUBPLOT_SIZE,SUBPLOT_SIZE), hadj=0)

subplot(fun = {
  with(mtreedat$data[which(mtreedat$data$genus %in% 'Sericomyia'),], plot(Comp.2 ~ Comp.1, pch=16, cex=0.5, col='red', tck=F,xlab='', ylab='', xaxt='n', yaxt='n', ylim=range(mtreedat$data$Comp.2[mtreedat$data$genus%in%genus6]), xlim=range(mtreedat$data$Comp.1[mtreedat$data$genus%in%genus6])));
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white");
  with(mtreedat$data[which(mtreedat$data$genus %in% 'Sericomyia'),], points(Comp.2 ~ Comp.1, pch=16, cex=0.5, col='red'))
}, x=R6, y=ySericomyia, size=c(SUBPLOT_SIZE,SUBPLOT_SIZE), hadj=0)

subplot(fun = {
  with(mtreedat$data[which(mtreedat$data$genus %in% 'Sphecomyia'),], plot(Comp.2 ~ Comp.1, pch=16, cex=0.5, col='red', tck=F,xlab='', ylab='', xaxt='n', yaxt='n', ylim=range(mtreedat$data$Comp.2[mtreedat$data$genus%in%genus6]), xlim=range(mtreedat$data$Comp.1[mtreedat$data$genus%in%genus6])));
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white");
  with(mtreedat$data[which(mtreedat$data$genus %in% 'Sphecomyia'),], points(Comp.2 ~ Comp.1, pch=16, cex=0.5, col='red'))
}, x=R6, y=ySphecomyia, size=c(SUBPLOT_SIZE,SUBPLOT_SIZE), hadj=0)

dev.off()

