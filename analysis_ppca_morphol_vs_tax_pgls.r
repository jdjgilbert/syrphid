

#### Analysis of morphology WRT different taxonomic levels

mtreedat$data$genus <- factor(sub('_.*','',row.names(mtreedat$data)))
mtreedat$data$tri <- syrph_tax$tri[match(mtreedat$data$genus, syrph_tax$genus)]
mtreedat$data$subtri <- syrph_tax$subtri[match(mtreedat$data$genus, syrph_tax$genus)]
mtreedat$data$subf <- syrph_tax$subf[match(mtreedat$data$tri, syrph_tax$tri)]

head(mtreedat$data)
#       SEX   WL   WW   HW   PL   FU  THW  THL  THH THVOL   T2   T3   T4  HTL LABR   PR   LL PS       Comp.1
#       Microdon_manitobensis/cot   M 6.67 2.58 3.38 2.05 0.66 2.60 3.71 2.88 27.78 4.40 4.40 3.85 2.30 0.51 0.42 0.76 33  0.036381238
#       Microdon_cothurnatus        F 6.68 2.63 3.41 2.09 0.70 2.40 3.41 2.64 21.61 4.22 4.12 3.80 2.06 0.47 0.45 0.84 31  0.069510146
#       Microdon_tristis            F 8.16 3.16 3.82 2.25 0.75 2.90 3.84 3.22 35.86 4.58 4.48 4.27 2.21 0.52 0.49 1.05 35 -0.090420434
#       Microdon_piperi             F 8.80 3.32 4.30 2.50 0.79 3.46 4.60 3.50 55.71 5.28 5.20 5.10 2.83 0.66 0.48 1.10 39 -0.260246527
#       Pipiza_crassipes            F 6.13 2.37 2.15 1.67 0.55 1.74 2.53 1.80  7.92 2.50 2.43 2.16 1.47 0.39 0.36 0.65 27  0.363129765
#       Pipiza_oregona              M 8.71 2.97 2.98 2.40 0.75 2.74 3.70 2.70 27.37 3.60 3.55 3.27 2.35 0.61 0.54 0.95 39  0.007285022
#       Comp.2     Comp.3        Comp.4    genus          tri       subtri          subf
#       Microdon_manitobensis/cot 0.2046104 0.17247761 -0.0681726194 Microdon Microdontini Microdontini Microdontinae
#       Microdon_cothurnatus      0.1817005 0.17604318 -0.0459042675 Microdon Microdontini Microdontini Microdontinae
#       Microdon_tristis          0.2322602 0.16757072 -0.0481638110 Microdon Microdontini Microdontini Microdontinae
#       Microdon_piperi           0.2773577 0.19277976 -0.0812345321 Microdon Microdontini Microdontini Microdontinae
#       Pipiza_crassipes          0.1265942 0.05213830 -0.0197206590   Pipiza     Pipizini     Pipizini     Syrphinae
#       Pipiza_oregona            0.1963788 0.04537503  0.0006184223   Pipiza     Pipizini     Pipizini     Syrphinae


pgls.subf <- pgls(Comp.2 ~ Comp.1 * subf + SEX - 1, data=mtreedat, lambda='ML')   	 ## 2014-09-04 
pgls.subf1 <- pgls(Comp.2 ~ Comp.1 + subf + SEX - 1, data=mtreedat, lambda=pgls.subf$param[2])
anova(pgls.subf, pgls.subf1)

##### OLD RESULTS AS OF 2014-09-04
#	Analysis of Variance Table
#	pgls: lambda = 0.96, delta = 1.00, kappa = 1.00
#	
#	Model 1: Comp.2 ~ Comp.1 * subf + SEX - 1
#	Model 2: Comp.2 ~ Comp.1 + subf + SEX - 1
#	  Res.Df    RSS Df Sum of Sq     F   Pr(>F)   
#	1    322 29.090                               
#	2    323 29.978 -1  -0.88814 9.831 0.001874 **
#	---
#	Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##### RESULTS AS OF 2015-12-22 - NEED TO DEBUG THIS IF NECESSARY
##   Analysis of Variance Table
##   pgls: lambda = 0.96, delta = 1.00, kappa = 1.00
##   
##   Model 1: Comp.2 ~ Comp.1 * subf + SEX - 1
##   Model 2: Comp.2 ~ Comp.1 + subf + SEX - 1
##   Res.Df    RSS Df Sum of Sq     F    Pr(>F)    
##   1    320 27.888                                 
##   2    322 29.673 -2   -1.7848 10.24 4.891e-05 ***

lrt(pgls.subf, pgls.subf1)
### 2014 result
# $LR
# [1] 9.834297
# 
# $DF
# [1] 1
# 
# $P
# [1] 0.001712877

### 2015/6 result
## $LR
## [1] 20.28536
## 
## $DF
## [1] 2
## 
## $P
## [1] 3.936315e-05
 
########### 2014-11-14: Head width as covariate
pgls.subf <- pgls(Comp.2 ~ HW * subf + SEX - 1, data=mtreedat, lambda='ML')   	 ## 2014-11-14
pgls.subf1 <- pgls(Comp.2 ~ HW + subf + SEX - 1, data=mtreedat, lambda=pgls.subf$param[2])
anova(pgls.subf, pgls.subf1)



#### Genus-level analysis #### 

genusdat <- subset(mtreedat$data, genus %in% names(a<-table(genus))[which(a>5)])
genusdat$genus <- factor(genusdat$genus)
genusdat$tri <- factor(genusdat$tri)
genusdat$binomial <- row.names(genusdat)
mtreedat.genus <- comparative.data(drop.tip(mtreedat$phy, tip=which(!mtreedat$phy$tip.label %in% row.names(genusdat))), genusdat, 'binomial')

pgls.genus <- pgls(Comp.2 ~ Comp.1 * genus + SEX - 1, data=mtreedat.genus, lambda='ML')   	 ## 2014-09-04 
pgls.genus1 <- pgls(Comp.2 ~ Comp.1 + genus + SEX - 1, data=mtreedat.genus, lambda=pgls.genus$param[2])
anova(pgls.genus, pgls.genus1)
# Analysis of Variance Table
# pgls: lambda = 0.00, delta = 1.00, kappa = 1.00
# 
# Model 1: Comp.2 ~ Comp.1 * genus + SEX - 1
# Model 2: Comp.2 ~ Comp.1 + genus + SEX - 1
# Res.Df     RSS  Df Sum of Sq      F    Pr(>F)    
# 1    156 0.29897                                   
# 2    174 0.49749 -18  -0.19852 5.7549 2.252e-10 ***

lrt(pgls.genus,pgls.genus1)
# $LR
# [1] 99.30155
# 
# $DF
# [1] 18
# 
# $P
# [1] 2.973177e-13

## EFfect of sex (2015-03-30)
pgls.genus2 <- pgls(Comp.2 ~ Comp.1 * genus - 1, data=mtreedat.genus, lambda='ML')   	 ## 2014-09-04 
lrt(pgls.genus,pgls.genus2)
# $LR
# [1] 8.839077
# 
# $DF
# [1] 1
# 
# $P
# [1] 0.002948482


########### 2014-11-14: Head width as covariate
pgls.genus <- pgls(Comp.2 ~ HW * genus + SEX - 1, data=mtreedat.genus, lambda='ML')   	 
pgls.genus1 <- pgls(Comp.2 ~ HW + genus + SEX - 1, data=mtreedat.genus, lambda=pgls.genus$param[2])
anova(pgls.genus, pgls.genus1)
#	Analysis of Variance Table
#	pgls: lambda = 0.00, delta = 1.00, kappa = 1.00
#	
#	Model 1: Comp.2 ~ HW * genus + SEX - 1
#	Model 2: Comp.2 ~ HW + genus + SEX - 1
#	  Res.Df     RSS  Df Sum of Sq      F    Pr(>F)    
#	1    156 0.34221                                   
#	2    174 0.50912 -18  -0.16691 4.2272 3.091e-07 ***
#	
lrt(pgls.genus,pgls.genus1)
#	$LR
#	[1] 77.46669
#	
#	$DF
#	[1] 18
#	
#	$P
#	[1] 2.368957e-09
#	

#### Subtri-level analysis #### 
 
subtridat <- subset(mtreedat$data, subtri %in% names(a<-table(mtreedat$data$subtri))[which(a>5)]) 	 ## 2014-09-04 
subtridat$subtri <- factor(subtridat$subtri)
subtridat$genus <- factor(subtridat$genus)
subtridat$binomial <- row.names(subtridat)
mtreedat.subtri <- comparative.data(drop.tip(mtreedat$phy, tip=which(!mtreedat$phy$tip.label %in% row.names(subtridat))), subtridat, 'binomial')

pgls.subtri <- pgls(Comp.2 ~ Comp.1 * subtri + SEX - 1, data=mtreedat.subtri, lambda='ML')  ## + SEX (if SEX important)
pgls.subtri1 <- pgls(Comp.2 ~ Comp.1 + subtri + SEX - 1, data=mtreedat.subtri, lambda=pgls.subtri$param[2])
pgls.subtri2 <- pgls(Comp.2 ~ Comp.1 * subtri - 1, data=mtreedat.subtri, lambda=pgls.subtri$param[2])
anova(pgls.subtri,pgls.subtri1)

#       Analysis of Variance Table
#       pgls: lambda = 0.73, delta = 1.00, kappa = 1.00
#       
#       Model 1: Comp.2 ~ Comp.1 * subtri + SEX
#       Model 2: Comp.2 ~ Comp.1 + subtri + SEX
#       Res.Df    RSS  Df Sum of Sq     F    Pr(>F)    
#       1    274 3.9072                                  
#       2    288 4.8451 -14  -0.93788 4.698 1.077e-07 ***
#         ---
#         Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

plot(resid(pgls.subtri)~fitted(pgls.subtri))
qqnorm(resid(pgls.subtri))  ## both lovely

lrt(pgls.subtri, pgls.subtri1)
#	$LR
#	[1] 65.61928
#	
#	$DF
#	[1] 14
#	
#	$P
#	[1] 1.185816e-08

lrt(pgls.subtri, pgls.subtri2)
#	$LR
#	[1] 9.271532
#	
#	$DF
#	[1] 1
#	
#	$P
#	[1] 0.00232743

########### 2014-11-14: Head width as covariate

pgls.subtri <- pgls(Comp.2 ~ HW * subtri + SEX - 1, data=mtreedat.subtri, lambda='ML')  ## + SEX (if SEX important)
pgls.subtri1 <- pgls(Comp.2 ~ HW + subtri + SEX - 1, data=mtreedat.subtri, lambda=pgls.subtri$param[2])
pgls.subtri2 <- pgls(Comp.2 ~ HW * subtri - 1, data=mtreedat.subtri, lambda=pgls.subtri$param[2])
anova(pgls.subtri,pgls.subtri1)
#	Analysis of Variance Table
#	pgls: lambda = 0.92, delta = 1.00, kappa = 1.00
#	
#	Model 1: Comp.2 ~ HW * subtri + SEX
#	Model 2: Comp.2 ~ HW + subtri + SEX
#	  Res.Df    RSS  Df Sum of Sq      F    Pr(>F)    
#	1    274 10.361                                   
#	2    288 12.630 -14   -2.2695 4.2872 7.286e-07 ***
lrt(pgls.subtri,pgls.subtri1)
#	$LR
#	[1] 60.41306
#	
#	$DF
#	[1] 14
#	
#	$P
#	[1] 9.928357e-08

lrt(pgls.subtri,pgls.subtri2)
#	$LR
#	[1] 9.479289
#	
#	$DF
#	[1] 1
#	
#	$P
#	[1] 0.002078045
  
#### Tribe-level analysis #### 

 tridat <- subset(mtreedat$data, tri %in% names(a<-table(mtreedat$data$tri))[which(a>5)])
tridat$tri <- factor(tridat$tri)
tridat$genus <- factor(tridat$genus)
tridat$binomial <- row.names(tridat)
mtreedat.tri <- comparative.data(drop.tip(mtreedat$phy, tip=which(!mtreedat$phy$tip.label %in% row.names(tridat))), tridat, 'binomial')

pgls.tri <- pgls(Comp.2 ~ Comp.1 * tri + SEX - 1, data=mtreedat.tri, lambda='ML') 	 ## 2014-09-04 
pgls.tri1 <- pgls(Comp.2 ~ Comp.1 + tri + SEX - 1, data=mtreedat.tri, lambda=pgls.tri$param[2])
pgls.tri2 <- pgls(Comp.2 ~ Comp.1 * tri - 1, data=mtreedat.tri, lambda=pgls.tri$param[2])
anova(pgls.tri, pgls.tri1)

##    Analysis of Variance Table
##    pgls: lambda = 0.93, delta = 1.00, kappa = 1.00
##    
##    Model 1: Comp.2 ~ Comp.1 * tri + SEX
##    Model 2: Comp.2 ~ Comp.1 + tri + SEX
##    Res.Df    RSS Df Sum of Sq      F    Pr(>F)    
##    1    295 18.597                                  
##    2    303 20.834 -8   -2.2364 4.4344 4.351e-05 ***
##      ---
##      Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
### 

lrt(pgls.tri, pgls.tri1)
#	$LR
#	[1] 35.65638
#	
#	$DF
#	[1] 8
#	
#	$P
#	[1] 2.029453e-05

lrt(pgls.tri, pgls.tri2)
#	$LR
#	[1] 5.884304
#	
#	$DF
#	[1] 1
#	
#	$P
#	[1] 0.01527643


########### 2014-11-14: Head width as covariate

pgls.tri <- pgls(Comp.2 ~ HW * tri + SEX - 1, data=mtreedat.tri, lambda='ML') 	 ## 2014-09-04 
pgls.tri1 <- pgls(Comp.2 ~ HW + tri + SEX - 1, data=mtreedat.tri, lambda=pgls.tri$param[2])
pgls.tri2 <- pgls(Comp.2 ~ HW * tri - 1, data=mtreedat.tri, lambda=pgls.tri$param[2])
anova(pgls.tri, pgls.tri1)
#	Analysis of Variance Table
#	pgls: lambda = 0.94, delta = 1.00, kappa = 1.00
#	
#	Model 1: Comp.2 ~ HW * tri + SEX
#	Model 2: Comp.2 ~ HW + tri + SEX
#	  Res.Df    RSS Df Sum of Sq      F   Pr(>F)   
#	1    295 21.988                                
#	2    303 23.732 -8   -1.7438 2.9244 0.003714 **

lrt(pgls.tri,pgls.tri1)
#	$LR
#	[1] 23.96384
#	
#	$DF
#	[1] 8
#	
#	$P
#	[1] 0.002324005

lrt(pgls.tri,pgls.tri2)
#	$LR
#	[1] 5.891827
#	
#	$DF
#	[1] 1
#	
#	$P
#	[1] 0.01521131

  
mtreedat.genus$data$subtri <- factor(mtreedat.genus$data$subtri)
mtreedat.genus$data$tri <- factor(mtreedat.genus$data$tri)


### Combine genus, tri and subtri into AIC analysis
pgls.genus <- pgls(Comp.2 ~ Comp.1 + SEX + genus + Comp.1:genus, data=mtreedat.genus, lambda='ML') 	 ## 2014-09-08 
pgls.genus1 <- pgls(Comp.2 ~ Comp.1 + SEX + genus, data=mtreedat.genus, lambda=pgls.genus$param[2]) 	 ## 2014-09-08 
pgls.tri <- pgls(Comp.2 ~ Comp.1 + SEX + tri + Comp.1:tri, data=mtreedat.genus, lambda='ML') 	 ## 2014-09-08 
pgls.tri1 <- pgls(Comp.2 ~ Comp.1 + SEX + tri, data=mtreedat.genus, lambda=pgls.tri$param[2]) 	 ## 2014-09-08 
pgls.subtri <- pgls(Comp.2 ~ Comp.1 + SEX + subtri + Comp.1:subtri, data=mtreedat.genus, lambda='ML') 	 ## 2014-09-08 
pgls.subtri1 <- pgls(Comp.2 ~ Comp.1 + SEX + subtri, data=mtreedat.genus, lambda=pgls.subtri$param[2]) 	 ## 2014-09-08 
pgls.subf <- pgls(Comp.2 ~ Comp.1 + SEX + subf + Comp.1:subf, data=mtreedat.genus, lambda='ML') 	 ## 2014-09-08 
pgls.subf1 <- pgls(Comp.2 ~ Comp.1 + SEX + subf, data=mtreedat.genus, lambda=pgls.subf$param[2]) 	 ## 2014-09-08 
pgls.null <- pgls(Comp.2 ~ Comp.1 + SEX, data=mtreedat.genus, lambda='ML') 	 ## 2014-09-08 

### Analyse again using HW as size measure
pgls.genus <- pgls(Comp.2 ~ HW + SEX + genus + HW:genus, data=mtreedat.genus, lambda='ML') 	 ## 2014-11-14 
pgls.genus1 <- pgls(Comp.2 ~ HW + SEX + genus, data=mtreedat.genus, lambda=pgls.genus$param[2]) 	 ## 2014-11-14 
pgls.tri <- pgls(Comp.2 ~ HW + SEX + tri + HW:tri, data=mtreedat.genus, lambda='ML') 	 ## 2014-11-14 
pgls.tri1 <- pgls(Comp.2 ~ HW + SEX + tri, data=mtreedat.genus, lambda=pgls.tri$param[2]) 	 ## 2014-11-14 
pgls.subtri <- pgls(Comp.2 ~ HW + SEX + subtri + HW:subtri, data=mtreedat.genus, lambda='ML') 	 ## 2014-11-14 
pgls.subtri1 <- pgls(Comp.2 ~ HW + SEX + subtri, data=mtreedat.genus, lambda=pgls.subtri$param[2]) 	 ## 2014-11-14 
pgls.subf <- pgls(Comp.2 ~ HW + SEX + subf + HW:subf, data=mtreedat.genus, lambda='ML') 	 ## 2014-11-14 
pgls.subf1 <- pgls(Comp.2 ~ HW + SEX + subf, data=mtreedat.genus, lambda=pgls.subf$param[2]) 	 ## 2014-11-14 
pgls.null <- pgls(Comp.2 ~ HW + SEX, data=mtreedat.genus, lambda='ML') 	 ## 2014-11-14 

### To paste into table
b<-(a <- AIC(pgls.genus,pgls.genus1,pgls.tri,pgls.tri1,pgls.subtri,pgls.subtri1,pgls.null))[order(a$AIC, decreasing=F),]  # pgls.subf,pgls.subf1,
b$daic <- b$AIC-min(b$AIC)
b$rla <- exp(-0.5 * b$daic)
b$wi <- b$rla / sum (b$rla)
b <- round(b, digits=4)

b
##    df       AIC    daic rla wi
##    pgls.genus   39 -613.7538  0.0000   1  1
##    pgls.subtri  25 -578.8006 34.9532   0  0
##    pgls.genus1  21 -572.2871 41.4667   0  0
##    pgls.null     3 -547.2225 66.5313   0  0
##    pgls.tri     15 -543.2291 70.5247   0  0
##    pgls.subtri1 14 -540.6730 73.0808   0  0
##    pgls.tri1     9 -537.6770 76.0768   0  0

### To paste into table
b<-(a <- AIC(pgls.genus,pgls.genus1,pgls.tri,pgls.tri1,pgls.subtri,pgls.subtri1,pgls.subf,pgls.subf1,pgls.null))[order(a$AIC, decreasing=F),]
b$daic <- b$AIC-min(b$AIC)
b$rla <- exp(-0.5 * b$daic)
b$wi <- b$rla / sum (b$rla)
b <- round(b, digits=4)


