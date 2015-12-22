
#### Karyotype-group analysis (PGLS and D-PGLS)


# Plot all genera with 4, 5, and 6 chromosomes

karnd$genus <- sub('_.*','',karnd$binomial)
nrow(a <- with(karnd, ftable(genus,Npair)))
head(a <- matrix(with(karnd, ftable(genus,sub('4.5','5', sub('5.5','6', sub('4.66666666666667','5', Npair))))), c(106,6)))
[,1] [,2] [,3] [,4] [,5]   ### 2014-12-01
[1,]    0    0    4    0    0
[2,]    0    7    7    0    0
[3,]    0    0    2    0    0
[4,]    0    0    0    0    1
[5,]    0    2    1    0    0
[6,]    0    0    4    0    0

row.names(a)<-sort(unique(karnd$genus))
colnames(a)<-sort(unique(sub('4.5','5', sub('5.5','6', sub('4.66666666666667','5', karnd$Npair)))))
a <- as.data.frame(a)

a$main <- apply(subset(a, select=c(1:5)), 1, function(x) names(a)[which.max(x)])
a$prop <- apply(subset(a, select=c(1:5)), 1, function(x) max(x)/sum(x))
a$occ <- as.vector(table(karnd$genus))

genus3 <- row.names(subset(a, main==3 & prop>0.60))
genus4 <- row.names(subset(a, main==4 & prop>0.60))
genus5 <- row.names(subset(a, main==5 & prop>0.60))
genus6 <- row.names(subset(a, main==6 & prop>0.60))
genus7 <- row.names(subset(a, main==7 & prop>0.60))

genusn <- table(mtreedat$data$genus)
genus4 <- genus4[which(genus4 %in% names(genusn[genusn>4]))]
genus5 <- genus5[which(genus5 %in% names(genusn[genusn>4]))]
genus6 <- genus6[which(genus6 %in% names(genusn[genusn>4]))]

mtreedat$data$gengrp <- 'AMB'
mtreedat$data$gengrp[mtreedat$data$genus %in% genus4] <- '4CH'
mtreedat$data$gengrp[mtreedat$data$genus %in% genus5] <- '5CH'
mtreedat$data$gengrp[mtreedat$data$genus %in% genus6] <- '6CH'

genusAMB <- unique(as.character(mtreedat$data$genus[which(mtreedat$data$gengrp=='AMB')]))

#with(subset(mtreedat$data, genus %in% genus3 & genus %in% names(genusn[genusn>3])), xyplot(Comp.2 ~ Comp.1 | genus)) ## only Baccha
#with(subset(mtreedat$data, genus %in% genus4 & genus %in% names(genusn[genusn>3])), xyplot(Comp.2 ~ Comp.1 | genus))
#with(subset(mtreedat$data, genus %in% genus5 & genus %in% names(genusn[genusn>3])), xyplot(Comp.2 ~ Comp.1 | genus))
#with(subset(mtreedat$data, genus %in% genus6 & genus %in% names(genusn[genusn>3])), xyplot(Comp.2 ~ Comp.1 | genus))
#with(subset(mtreedat$data, genus %in% genus7 & genus %in% names(genusn[genusn>3])), xyplot(Comp.2 ~ Comp.1 | genus))



### is belonging to the genus-karyotype group associated with allometry?  YES.

Y <- as.matrix(mtreedat$data[,2:grep('PS', names(mtreedat$data)),])  ### Y is a matrix of all morphological data

procD.pgls(Y ~ mtreedat$data$SEX + mtreedat$data$Comp.1 * mtreedat$data$gengrp, phy=mtreedat$phy, iter=1000)
#	                                           df       SS       MS     Rsq        F   P.val
#	mtreedat$data$SEX                           1    96119    96119 0.00373   4.8433 0.21678
#	mtreedat$data$Comp.1                        1 18847893 18847893 0.73185 949.7215 0.00100
#	mtreedat$data$gengrp                        3     8337     2779 0.00032   0.1400 0.49251
#	mtreedat$data$Comp.1:mtreedat$data$gengrp   3   490639   163546 0.01905   8.2409 0.00999
#	Residuals                                 318  6310934    19846     


### now look only among gengrps (ie. remove the AMBIGUOUS group)

nrow(mtreedat.gengrp.data <- subset(mtreedat$data, genus %in% c(genus4, genus5, genus6) & genus %in% names(genusn[genusn>4])))
# 252
mtreedat.gengrp.data$binomial <- row.names(mtreedat.gengrp.data)
mtreedat.gengrp.data$gengrp <- '4CH'
mtreedat.gengrp.data$gengrp[mtreedat.gengrp.data$genus %in% genus5] <- '5CH'
mtreedat.gengrp.data$gengrp[mtreedat.gengrp.data$genus %in% genus6] <- '6CH'
mtreedat.gengrp <- comparative.data(phy=drop.tip(mtreedat$phy, tip=which(!mtreedat$phy$tip.label %in% row.names(mtreedat.gengrp.data))),data=mtreedat.gengrp.data,'binomial')
mtreedat.gengrp$data$genus <- factor(mtreedat.gengrp$data$genus)
Y <- as.matrix(mtreedat.gengrp$data[,2:grep('PS', names(mtreedat.gengrp$data)),])  ### Y is a matrix of all morphological data


##########
procD.pgls(Y ~ mtreedat.gengrp$data$SEX + mtreedat.gengrp$data$Comp.1 * mtreedat.gengrp$data$gengrp, phy=mtreedat.gengrp$phy, iter=1000)
#	                                                         df       SS       MS     Rsq        F    P.val
#	mtreedat.gengrp$data$SEX                                  1   116570   116570 0.00834   8.9981 0.086913
#	mtreedat.gengrp$data$Comp.1                               1 10853898 10853898 0.77611 837.8178 0.000999
#	mtreedat.gengrp$data$gengrp                               2     2223     1112 0.00016   0.0858 0.025974
#	mtreedat.gengrp$data$Comp.1:mtreedat.gengrp$data$gengrp   2   447167   223584 0.03197  17.2585 0.004995
#	Residuals                                               198  2565082    12955 

procD.pgls(Y[,-3] ~ mtreedat.gengrp$data$SEX + mtreedat.gengrp$data$HW * mtreedat.gengrp$data$gengrp, phy=mtreedat.gengrp$phy, iter=1000)
#	                                                     df       SS       MS     Rsq        F    P.val
#	mtreedat.gengrp$data$SEX                              1   116445   116445 0.00833   6.7698 0.110889
#	mtreedat.gengrp$data$HW                               1 10140456 10140456 0.72539 589.5381 0.000999
#	mtreedat.gengrp$data$gengrp                           2     2162     1081 0.00015   0.0628 0.058941
#	mtreedat.gengrp$data$HW:mtreedat.gengrp$data$gengrp   2   314436   157218 0.02249   9.1402 0.009990
#	Residuals                                           198  3405735    17201 
##########



mtreedat$data$gengrp <- factor(mtreedat$data$gengrp, levels=c('AMB','4CH','5CH','6CH'))

contrasts(mtreedat$data$gengrp) <- cbind(c(-3,1,1,1), c(0, -1, -1, 2), c(0, -1, 1, 0))
## 1: NONE vs rest
## 2: (4,5) vs 6
## 3: 4 vs 5

gls1 <- gls(Comp.2 ~ SEX + Comp.1 * gengrp, data=mtreedat$data, correlation=corPagel(1, mtreedat$phy),  method='ML')
gls2 <- gls(Comp.2 ~ SEX + Comp.1 + gengrp, data=mtreedat$data, correlation=corPagel(1, mtreedat$phy),  method='ML')
anova(gls1, gls2)
#	     Model df       AIC       BIC   logLik   Test  L.Ratio p-value
#	gls1     1 11 -651.8933 -610.2038 336.9467                        
#	gls2     2  8 -643.1007 -612.7810 329.5503 1 vs 2 14.79265   0.002

summary(gls1)
#	Generalized least squares fit by maximum likelihood
#	  Model: Comp.2 ~ SEX + Comp.1 * gengrp 
#	  Data: mtreedat$data 
#	        AIC       BIC   logLik
#	  -651.8933 -610.2038 336.9467
#	
#	Correlation Structure: corPagel
#	 Formula: ~1 
#	 Parameter estimate(s):
#	   lambda 
#	0.9479143 
#	
#	Coefficients:
#	                     Value  Std.Error    t-value p-value
#	(Intercept)     0.05841204 0.13545098  0.4312412  0.6666
#	SEXM           -0.02768937 0.00952568 -2.9068141  0.0039
#	Comp.1         -0.02317998 0.02333561 -0.9933307  0.3213
#	gengrp1         0.00247768 0.00453752  0.5460416  0.5854
#	gengrp2         0.00173240 0.01115484  0.1553051  0.8767
#	gengrp3        -0.00542963 0.02083100 -0.2606517  0.7945
#	Comp.1:gengrp1  0.00820378 0.00872058  0.9407373  0.3476
#	Comp.1:gengrp2  0.04489863 0.01588331  2.8267798  0.0050
#	Comp.1:gengrp3 -0.01507999 0.03923345 -0.3843656  0.7010
#	
#	 Correlation: 
#	               (Intr) SEXM   Comp.1 gngrp1 gngrp2 gngrp3 Cm.1:1 Cm.1:2
#	SEXM           -0.020                                                 
#	Comp.1         -0.011 -0.015                                          
#	gengrp1         0.070 -0.058 -0.057                                   
#	gengrp2        -0.012 -0.035  0.051 -0.171                            
#	gengrp3         0.025 -0.037  0.057  0.270 -0.236                     
#	Comp.1:gengrp1  0.001  0.035  0.447  0.049  0.047  0.129              
#	Comp.1:gengrp2  0.013 -0.107 -0.493  0.081  0.125 -0.069 -0.474       
#	Comp.1:gengrp3  0.012  0.032  0.190  0.108 -0.037  0.018  0.278 -0.260
#	
#	Standardized residuals:
#	       Min         Q1        Med         Q3        Max 
#	-2.8406289 -0.4825859 -0.1563412  0.2328058  0.8126356 
#	
#	Residual standard error: 0.2765485 
#	Degrees of freedom: 327 total; 318 residual

#### USING HW AS SIZE MEASURE 2014-11-14
gls1 <- gls(Comp.2 ~ SEX + HW * gengrp, data=mtreedat$data, correlation=corPagel(1, mtreedat$phy),  method='ML')
gls2 <- gls(Comp.2 ~ SEX + HW + gengrp, data=mtreedat$data, correlation=corPagel(1, mtreedat$phy),  method='ML')
anova(gls1, gls2)
#	     Model df       AIC       BIC   logLik   Test  L.Ratio p-value
#	gls1     1 11 -658.1270 -616.4375 340.0635                        
#	gls2     2  8 -648.4915 -618.1718 332.2457 1 vs 2 15.63555  0.0013

summary(gls1)
#	Generalized least squares fit by maximum likelihood
#	  Model: Comp.2 ~ SEX + HW * gengrp 
#	  Data: mtreedat$data 
#	       AIC       BIC   logLik
#	  -658.127 -616.4375 340.0635
#	
#	Correlation Structure: corPagel
#	 Formula: ~1 
#	 Parameter estimate(s):
#	   lambda 
#	0.9576604 
#	
#	Coefficients:
#	                  Value  Std.Error    t-value p-value
#	(Intercept) -0.00620554 0.14730658 -0.0421267  0.9664
#	SEXM        -0.02528181 0.00927672 -2.7252963  0.0068
#	HW           0.02130025 0.00973029  2.1890665  0.0293
#	gengrp1      0.01951282 0.01138062  1.7145653  0.0874
#	gengrp2      0.04691488 0.02363625  1.9848702  0.0480
#	gengrp3     -0.00875540 0.05247380 -0.1668527  0.8676
#	HW:gengrp1  -0.00538309 0.00345309 -1.5589187  0.1200
#	HW:gengrp2  -0.01599575 0.00665569 -2.4033188  0.0168
#	HW:gengrp3   0.00087716 0.01625469  0.0539633  0.9570
#	
#	 Correlation: 
#	           (Intr) SEXM   HW     gngrp1 gngrp2 gngrp3 HW:gn1 HW:gn2
#	SEXM       -0.021                                                 
#	HW         -0.185  0.011                                          
#	gengrp1     0.102  0.001 -0.417                                   
#	gengrp2    -0.066 -0.076  0.367 -0.351                            
#	gengrp3     0.048 -0.036 -0.159  0.296 -0.226                     
#	HW:gengrp1 -0.092 -0.027  0.506 -0.913  0.388 -0.247              
#	HW:gengrp2  0.081  0.079 -0.485  0.402 -0.877  0.206 -0.493       
#	HW:gengrp3 -0.038  0.020  0.156 -0.243  0.180 -0.913  0.236 -0.206
#	
#	Standardized residuals:
#	       Min         Q1        Med         Q3        Max 
#	-2.6023649 -0.4580260 -0.1054674  0.2172115  0.7437722 
#	
#	Residual standard error: 0.2942698 
#	Degrees of freedom: 327 total; 318 residual


##### REPEAT WITHIN ERISTALINAE - OBVIOUSLY INNOVATION OF 6 CHROMOSOMES WITHIN THIS SUBF HAS CHANGED ALLOMETRY
eris.gengrp.phy1 <- drop.tip(mtreedat$phy, tip=row.names(subset(mtreedat$data, subf!='Eristalinae')))
eris.gengrp.dat1 <- within((a <- subset(mtreedat$data, subf=='Eristalinae')), binomial <- row.names(a))
kmtd.gengrp.eris <- comparative.data(eris.gengrp.phy1, eris.gengrp.dat1, 'binomial')
kmtd.gengrp.eris$data$genus <- factor(kmtd.gengrp.eris$data$genus)


gls1 <- gls(Comp.2 ~ SEX + Comp.1 * gengrp, data=kmtd.gengrp.eris$data, correlation=corPagel(1, kmtd.gengrp.eris$phy),  method='ML')
gls2 <- gls(Comp.2 ~ SEX + Comp.1 + gengrp, data=kmtd.gengrp.eris$data, correlation=corPagel(1, kmtd.gengrp.eris$phy),  method='ML')
anova(gls1, gls2)
#	     Model df       AIC       BIC   logLik   Test  L.Ratio p-value
#	gls1     1 11 -375.4301 -339.2591 198.7150                        
#	gls2     2  8 -374.0110 -347.7048 195.0055 1 vs 2 7.419103  0.0597

summary(gls1)
#	Generalized least squares fit by maximum likelihood
#	  Model: Comp.2 ~ SEX + Comp.1 * gengrp 
#	  Data: kmtd.gengrp.eris$data 
#	        AIC       BIC   logLik
#	  -373.7699 -337.5989 197.8849
#	
#	Correlation Structure: corPagel
#	 Formula: ~1 
#	 Parameter estimate(s):
#	   lambda 
#	0.8543994 
#	
#	Coefficients:
#	                     Value  Std.Error    t-value p-value
#	(Intercept)    -0.04870069 0.08058528 -0.6043374  0.5463
#	SEXM           -0.02322057 0.01391756 -1.6684371  0.0969
Comp.1         -0.00365550 0.03462865 -0.1055629  0.9160
gengrp1         0.00858646 0.00762631  1.1258999  0.2616
gengrp2        -0.01143833 0.01529354 -0.7479189  0.4554
gengrp3         0.01490682 0.03907555  0.3814872  0.7033
Comp.1:gengrp1 -0.01178735 0.01277528 -0.9226684  0.3574
Comp.1:gengrp2  0.05192335 0.02431455  2.1354850  0.0340
Comp.1:gengrp3  0.03466790 0.06513753  0.5322262  0.5952

Correlation: 
  (Intr) SEXM   Comp.1 gngrp1 gngrp2 gngrp3 Cm.1:1 Cm.1:2
SEXM           -0.060                                                 
Comp.1          0.048 -0.006                                          
gengrp1         0.224 -0.110  0.154                                   
gengrp2        -0.131 -0.007 -0.120 -0.575                            
gengrp3        -0.128 -0.013 -0.270 -0.332  0.238                     
Comp.1:gengrp1  0.071 -0.013  0.675  0.333 -0.229 -0.339              
Comp.1:gengrp2 -0.054 -0.124 -0.672 -0.228  0.347  0.358 -0.700       
Comp.1:gengrp3 -0.101  0.177  0.028 -0.319  0.296  0.283  0.013 -0.043

Standardized residuals:
  Min         Q1        Med         Q3        Max 
-3.3560913 -0.1187527  0.1590451  0.6914117  1.8807331 

Residual standard error: 0.1835567 
Degrees of freedom: 198 total; 189 residual



###### FIRST TWO CHARACTERS OF EACH LINE HAVE BEEN REMOVED FROM HEREON DOWN - FIX THIS !!


# Using head width
s1 <- gls(Comp.2 ~ SEX + HW * gengrp, data=kmtd.gengrp.eris$data, correlation=corPagel(1, kmtd.gengrp.eris$phy),  method='ML')
s2 <- gls(Comp.2 ~ SEX + HW + gengrp, data=kmtd.gengrp.eris$data, correlation=corPagel(1, kmtd.gengrp.eris$phy),  method='ML')
ova(gls1, gls2)
Model df       AIC       BIC   logLik   Test  L.Ratio p-value
gls1     1 11 -373.4471 -337.2762 197.7236                        
gls2     2  8 -371.5116 -345.2055 193.7558 1 vs 2 7.935478  0.0474

mmary(gls1)
Generalized least squares fit by maximum likelihood
Model: Comp.2 ~ SEX + HW * gengrp 
Data: kmtd.gengrp.eris$data 
AIC       BIC   logLik
-373.4471 -337.2762 197.7236

Correlation Structure: corPagel
Formula: ~1 
Parameter estimate(s):
  lambda 
0.8562133 

Coefficients:
  Value  Std.Error    t-value p-value
(Intercept) -0.09264281 0.09342292 -0.9916497  0.3226
SEXM        -0.02108930 0.01367938 -1.5416852  0.1248
HW           0.01780985 0.01514137  1.1762372  0.2410
gengrp1      0.01236813 0.02008824  0.6156902  0.5388
gengrp2      0.03491966 0.03748049  0.9316755  0.3527
gengrp3      0.02405516 0.10050093  0.2393526  0.8111
HW:gengrp1   0.00053202 0.00548191  0.0970504  0.9228
HW:gengrp2  -0.01949737 0.01025206 -1.9018011  0.0587
HW:gengrp3  -0.00812211 0.02824761 -0.2875328  0.7740

Correlation: 
  (Intr) SEXM   HW     gngrp1 gngrp2 gngrp3 HW:gn1 HW:gn2
SEXM       -0.059                                                 
HW         -0.501  0.022                                          
gengrp1     0.390 -0.075 -0.572                                   
gengrp2    -0.389 -0.028  0.616 -0.728                            
gengrp3    -0.194  0.108  0.151 -0.339  0.312                     
HW:gengrp1 -0.378  0.038  0.668 -0.932  0.704  0.219              
HW:gengrp2  0.391  0.044 -0.721  0.683 -0.923 -0.198 -0.750       
HW:gengrp3  0.119 -0.134 -0.045  0.198 -0.183 -0.926 -0.074  0.054

Standardized residuals:
  Min         Q1        Med         Q3        Max 
-3.4208351 -0.1382392  0.2030174  0.6263538  1.5899620 

Residual standard error: 0.1845612 
Degrees of freedom: 198 total; 189 residual



<- as.matrix(kmtd.gengrp.eris$data[,2:grep('PS', names(kmtd.gengrp.eris$data)),])  ### Y is a matrix of all morphological data
ocD.pgls(Y ~ kmtd.gengrp.eris$data$SEX + kmtd.gengrp.eris$data$Comp.1 * kmtd.gengrp.eris$data$gengrp, phy=kmtd.gengrp.eris$phy, iter=1000)  ##
df       SS       MS     Rsq        F   P.val
kmtd.gengrp.eris$data$SEX                                   1   162570   162570 0.00713   6.7901 0.13586
kmtd.gengrp.eris$data$Comp.1                                1 17748306 17748306 0.77805 741.2910 0.00100
kmtd.gengrp.eris$data$gengrp                                3     8194     2731 0.00036   0.1141 0.64935
kmtd.gengrp.eris$data$Comp.1:kmtd.gengrp.eris$data$gengrp   3   367010   122337 0.01609   5.1096 0.02198
Residuals                                                 189  4525119    23942         

repeat using HW
ocD.pgls(Y[,-3] ~ kmtd.gengrp.eris$data$SEX + kmtd.gengrp.eris$data$HW * kmtd.gengrp.eris$data$gengrp, phy=kmtd.gengrp.eris$phy, iter=1000)  ##
df       SS       MS     Rsq        F    P.val
kmtd.gengrp.eris$data$SEX                               1   162469   162469 0.00712   5.9242 0.151848
kmtd.gengrp.eris$data$HW                                1 17084923 17084923 0.74920 622.9780 0.000999
kmtd.gengrp.eris$data$gengrp                            3   126593    42198 0.00555   1.5387 0.022977
kmtd.gengrp.eris$data$HW:kmtd.gengrp.eris$data$gengrp   3   246873    82291 0.01083   3.0006 0.112887
Residuals                                             189  5183249    27425       

## EXCLUDING AMBIGUOUS GENERA
ow(eris.gengrp.data <- subset(kmtd.gengrp.eris$data, gengrp!='AMB'))
132
is.gengrp.data$binomial <- row.names(eris.gengrp.data)
is.gengrp.data$gengrp <- '4CH'
is.gengrp.data$gengrp[eris.gengrp.data$genus %in% genus5] <- '5CH'
is.gengrp.data$gengrp[eris.gengrp.data$genus %in% genus6] <- '6CH'
is.gengrp <- comparative.data(phy=drop.tip(mtreedat$phy, tip=which(!mtreedat$phy$tip.label %in% row.names(eris.gengrp.data))),data=eris.gengrp.data,'binomial')
is.gengrp$data$genus <- factor(eris.gengrp$data$genus)
<- as.matrix(eris.gengrp$data[,2:grep('PS', names(eris.gengrp$data)),])  ### Y is a matrix of all morphological data


########
ocD.pgls(Y ~ eris.gengrp$data$SEX + eris.gengrp$data$Comp.1 * eris.gengrp$data$gengrp, phy=eris.gengrp$phy, iter=1000)
df       SS       MS     Rsq        F    P.val
eris.gengrp$data$SEX                              1    94195    94195 0.00733   6.1219 0.118881
eris.gengrp$data$Comp.1                           1 10609255 10609255 0.82528 689.5196 0.000999
eris.gengrp$data$gengrp                           2     2480     1240 0.00019   0.0806 0.033966
eris.gengrp$data$Comp.1:eris.gengrp$data$gengrp   2   226047   113023 0.01758   7.3457 0.015984
Residuals                                       125  1923305    15386     

using HW
ocD.pgls(Y[,-3] ~ eris.gengrp$data$SEX + eris.gengrp$data$HW * eris.gengrp$data$gengrp, phy=eris.gengrp$phy, iter=1000)


