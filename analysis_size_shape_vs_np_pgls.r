
## PGLS analysis of size-shape relationship WRT NP (or KPPCA1)

### Does np affect Comp.1/Comp.2 relationship? YES - whether weighted or not

# unweighted
gls1 <- gls(Comp.2 ~ SEX + Comp.1 * np, data=kmtreedat$data, correlation=corPagel(1, kmtreedat$phy), method='ML')  	 ## 2014-11-14
gls2 <- gls(Comp.2 ~ SEX + Comp.1 + np, data=kmtreedat$data, correlation=corPagel(1, kmtreedat$phy), method='ML')
anova(gls1,gls2)
#	     Model df       AIC       BIC   logLik   Test  L.Ratio p-value
#	gls1     1  7 -306.8930 -285.2798 160.4465                        
#	gls2     2  6 -299.7347 -281.2091 155.8674 1 vs 2 9.158271  0.0025
lrt(gls1,gls2)
#	$LR
#	[1] 9.158271
#	
#	$DF
#	[1] 1
#	
#	$P
#	[1] 0.002475964

summary(gls1) ### for Supplementary table
#	Generalized least squares fit by maximum likelihood
#	  Model: Comp.2 ~ SEX + Comp.1 * np 
#	  Data: kmtreedat$data 
#	       AIC       BIC   logLik
#	  -306.893 -285.2798 160.4465
#	
#	Correlation Structure: corPagel
#	 Formula: ~1 
#	 Parameter estimate(s):
#	   lambda 
#	0.9656051 
#	
#	Coefficients:
#	                 Value  Std.Error   t-value p-value
#	(Intercept)  0.2306645 0.15801543  1.459760  0.1464
#	SEXM        -0.0234078 0.01341633 -1.744721  0.0830
#	Comp.1      -0.3645055 0.12349113 -2.951674  0.0036
#	np          -0.0292319 0.01326114 -2.204326  0.0290
#	Comp.1:np    0.0684999 0.02264621  3.024784  0.0029
#	
#	 Correlation: 
#	          (Intr) SEXM   Comp.1 np    
#	SEXM       0.054                     
#	Comp.1     0.083  0.026              
#	np        -0.442 -0.158 -0.261       
#	Comp.1:np -0.080 -0.053 -0.977  0.259
#	
#	Standardized residuals:
#	        Min          Q1         Med          Q3         Max 
#	-2.81665829 -0.51238534 -0.12475744  0.08308796  0.98430846 
#	
#	Residual standard error: 0.2857265 
#	Degrees of freedom: 162 total; 157 residual

### Using HW for size
gls1 <- gls(Comp.2 ~ SEX + HW * np, data=kmtreedat$data, correlation=corPagel(1, kmtreedat$phy), method='ML')  	 ## 2014-11-14
gls2 <- gls(Comp.2 ~ SEX + HW + np, data=kmtreedat$data, correlation=corPagel(1, kmtreedat$phy), method='ML')
anova(gls1,gls2)
#	     Model df       AIC       BIC   logLik   Test  L.Ratio p-value
#	gls1     1  7 -308.6332 -287.0200 161.3166                        
#	gls2     2  6 -301.1185 -282.5929 156.5592 1 vs 2 9.514676   0.002
lrt(gls1,gls2)
#	$LR
#	[1] 9.514676
#	
#	$DF
#	[1] 1
#	
#	$P
#	[1] 0.00203835

summary(gls1) # for Supplementary table
#	Generalized least squares fit by maximum likelihood
#	  Model: Comp.2 ~ SEX + HW * np 
#	  Data: kmtreedat$data 
#	        AIC     BIC   logLik
#	  -308.6332 -287.02 161.3166
#	
#	Correlation Structure: corPagel
#	 Formula: ~1 
#	 Parameter estimate(s):
#	  lambda 
#	0.968117 
#	
#	Coefficients:
#	                  Value  Std.Error   t-value p-value
#	(Intercept) -0.26236411 0.22689031 -1.156348  0.2493
#	SEXM        -0.02251247 0.01323211 -1.701352  0.0909
#	HW           0.16823924 0.05146739  3.268851  0.0013
#	np           0.05885754 0.03394964  1.733672  0.0849
#	HW:np       -0.02959385 0.00957941 -3.089319  0.0024
#	
#	 Correlation: 
#	      (Intr) SEXM   HW     np    
#	SEXM   0.070                     
#	HW    -0.712 -0.050              
#	np    -0.757 -0.124  0.911       
#	HW:np  0.691  0.074 -0.977 -0.927
#	
#	Standardized residuals:
#	        Min          Q1         Med          Q3         Max 
#	-2.77865308 -0.55470132 -0.17233004  0.06573653  1.01686183 
#	
#	Residual standard error: 0.2901723 
#	Degrees of freedom: 162 total; 157 residual

### randomization test
res.gls1 <- res.gls2 <- list()
for (i in 1:100) {
  kmtreedat$data$npx <- KK[[i]]
  res.gls1[[i]] <- gls(Comp.2 ~ Comp.1 * npx, data=kmtreedat$data, correlation=corPagel(1, kmtreedat$phy), method='ML')	 ## 2014-09-09 
  res.gls2[[i]] <- gls(Comp.2 ~ Comp.1 + npx, data=kmtreedat$data, correlation=corPagel(1, kmtreedat$phy), method='ML')
  warning(i, immediate.=T)
}
res.gls.anova <- sapply(1:100, function(x) anova(res.gls1[[x]], res.gls2[[x]])$L.Ratio)
quantile(res.gls.anova[2,], c(0.05, 0.95))
5%       95% 
  0.0207576 5.1214590   ### 8.48 is well outside this

### np has no effect upon Comp.3
gls1 <- gls(Comp.3 ~ Comp.1 * np, data=kmtreedat$data, correlation=corPagel(1, kmtreedat$phy), method='ML')  	 ## 2014-09-09 
gls2 <- gls(Comp.3 ~ Comp.1 + np, data=kmtreedat$data, correlation=corPagel(1, kmtreedat$phy), method='ML')
anova(gls1,gls2)
#	     Model df       AIC       BIC   logLik   Test  L.Ratio p-value
#	gls1     1  6 -452.9277 -434.4021 232.4639                        
#	gls2     2  5 -454.8018 -439.3638 232.4009 1 vs 2 0.125958  0.7227
## or Comp.4
gls1 <- gls(Comp.4 ~ Comp.1 * np, data=kmtreedat$data, correlation=corPagel(1, kmtreedat$phy), method='ML')  	 ## 2014-09-09 
gls2 <- gls(Comp.4 ~ Comp.1 + np, data=kmtreedat$data, correlation=corPagel(1, kmtreedat$phy), method='ML')
anova(gls1,gls2)
#	     Model df       AIC       BIC   logLik   Test   L.Ratio p-value
#	gls1     1  6 -513.3998 -494.8742 262.6999                         
#	gls2     2  5 -514.7000 -499.2620 262.3500 1 vs 2 0.6997507  0.4029

### Look at KPPCA1



## Using PGLS ## 2014-09-02
#	pgls1 <- pgls(Comp.2 ~ Comp.1 * np, data=kmtreedat, lambda='ML')  ### Fails with weird message: Problem with optim:52ERROR: ABNORMAL_TERMINATION_IN_LNSRCH
#	pgls2 <- pgls(Comp.2 ~ Comp.1 + np, data=kmtreedat, lambda=pgls1$param[2])
#	anova(pgls1,pgls2)
#	lrt(pgls1,pgls2)

## Variance weighted by np

vf <- varFixed(~np)  ## look at variance, whether this increases with npair (it does)
vf <- Initialize(vf, kmtreedat$data)
vf1 <- varIdent(~np)  ## look at variance, whether this increases with npair (it does)
vf1 <- Initialize(vf1, kmtreedat$data)

gls1w <- gls(Comp.2 ~ Comp.1 * np, data=kmtreedat$data, correlation=corPagel(1, kmtreedat$phy), weights=vf, method='ML')
gls2w <- gls(Comp.2 ~ Comp.1 + np, data=kmtreedat$data, correlation=corPagel(1, kmtreedat$phy), weights=vf, method='ML')
anova(gls1w,gls2w)
#	      Model df       AIC       BIC   logLik   Test  L.Ratio p-value
#	gls1w     1  6 -302.1294 -283.6039 157.0647                        
#	gls2w     2  5 -296.6818 -281.2438 153.3409 1 vs 2 7.447636  0.0064

### Does weighting the variance matter? NO - unweighted model has the lowest AIC by >3 points
anova(gls1w,gls1)
Model df       AIC       BIC   logLik
gls1w     1  6 -302.1294 -283.6039 157.0647
gls1      2  6 -305.7825 -287.2569 158.8913

kmtreedat$data$genus <- factor(sub('_.*','',row.names(kmtreedat$data)))
kmtreedat$data$tri <- syrph_tax$tri[match(kmtreedat$data$genus, syrph_tax$genus)]
kmtreedat$data$subtri <- syrph_tax$subtri[match(kmtreedat$data$genus, syrph_tax$genus)]
kmtreedat$data$subf <- syrph_tax$subf[match(kmtreedat$data$tri, syrph_tax$tri)]




## Does this pattern change among subfamilies?  Include subf as a predictor
gls1 <- gls(Comp.2 ~ Comp.1 * np * subf, data=kmtreedat$data, correlation=corPagel(1, kmtreedat$phy), weights=vf, method='ML')
gls2 <- gls(Comp.2 ~ (Comp.1 + np + subf)^2, data=kmtreedat$data, correlation=corPagel(1, kmtreedat$phy), weights=vf, method='ML')
gls3 <- update(gls2, ~.-Comp.1:np)
gls4 <- update(gls2, ~.-Comp.1:subf)
gls5 <- update(gls2, ~.-np:subf)
anova(gls1,gls2)
#	     Model df       AIC       BIC   logLik   Test   L.Ratio p-value
#	gls1     1 10 -254.9979 -224.8915 137.4990                         
#	gls2     2  9 -256.6656 -229.5699 137.3328 1 vs 2 0.3323049  0.5643

anova(gls2,gls3)
#	     Model df       AIC       BIC   logLik   Test  L.Ratio p-value
#	gls2     1  9 -256.6656 -229.5699 137.3328                        
#	gls3     2  8 -252.6640 -228.5789 134.3320 1 vs 2 6.001602  0.0143
anova(gls2,gls4)
#	     Model df       AIC       BIC   logLik   Test  L.Ratio p-value
#	gls2     1  9 -256.6656 -229.5699 137.3328                        
#	gls4     2  8 -255.7958 -231.7107 135.8979 1 vs 2 2.869842  0.0903

anova(gls2,gls5)
#	     Model df       AIC       BIC   logLik   Test  L.Ratio p-value
#	gls2     1  9 -256.6656 -229.5699 137.3328                        
#	gls5     2  8 -257.6655 -233.5805 136.8328 1 vs 2 1.000058  0.3173

gls6 <- update(gls5, ~.-Comp.1:subf)
anova(gls5, gls6)
#	     Model df       AIC       BIC   logLik   Test  L.Ratio p-value
#	gls5     1  8 -257.6655 -233.5805 136.8328                        
#	gls6     2  7 -256.2829 -235.2085 135.1414 1 vs 2 3.382634  0.0659  .

gls7 <- update(gls6, ~.-subf)
#	     Model df       AIC       BIC   logLik   Test   L.Ratio p-value
#	gls6     1  7 -256.2829 -235.2085 135.1414                         
#	gls7     2  6 -257.9115 -239.8477 134.9558 1 vs 2 0.3713841  0.5423

## redo using AIC
gls1 <- gls(Comp.2 ~ Comp.1 * np * subf, data=kmtreedat$data, correlation=corPagel(1, kmtreedat$phy), weights=vf, method='ML')
gls2 <- gls(Comp.2 ~ Comp.1 + np + subf + Comp.1:np + Comp.1:subf + np:subf, data=kmtreedat$data, correlation=corPagel(1, kmtreedat$phy), weights=vf, method='ML')
gls3 <- gls(Comp.2 ~ Comp.1 + np + subf + Comp.1:np + Comp.1:subf          , data=kmtreedat$data, correlation=corPagel(1, kmtreedat$phy), weights=vf, method='ML')
gls4 <- gls(Comp.2 ~ Comp.1 + np + subf + Comp.1:np +               np:subf, data=kmtreedat$data, correlation=corPagel(1, kmtreedat$phy), weights=vf, method='ML')
gls5 <- gls(Comp.2 ~ Comp.1 + np + subf +             Comp.1:subf + np:subf, data=kmtreedat$data, correlation=corPagel(1, kmtreedat$phy), weights=vf, method='ML')
gls6 <- gls(Comp.2 ~ Comp.1 + np + subf + Comp.1:np                        , data=kmtreedat$data, correlation=corPagel(1, kmtreedat$phy), weights=vf, method='ML')
gls7 <- gls(Comp.2 ~ Comp.1 + np        + Comp.1:np                        , data=kmtreedat$data, correlation=corPagel(1, kmtreedat$phy), weights=vf, method='ML')
gls8 <- gls(Comp.2 ~ Comp.1 + np + subf             + Comp.1:subf          , data=kmtreedat$data, correlation=corPagel(1, kmtreedat$phy), weights=vf, method='ML')
gls9 <- gls(Comp.2 ~ Comp.1      + subf             + Comp.1:subf          , data=kmtreedat$data, correlation=corPagel(1, kmtreedat$phy), weights=vf, method='ML')
gls10<- gls(Comp.2 ~ Comp.1 + np + subf                           + np:subf, data=kmtreedat$data, correlation=corPagel(1, kmtreedat$phy), weights=vf, method='ML')
gls11<- gls(Comp.2 ~          np + subf                           + np:subf, data=kmtreedat$data, correlation=corPagel(1, kmtreedat$phy), weights=vf, method='ML')
gls12<- gls(Comp.2 ~ Comp.1 + np + subf                                    , data=kmtreedat$data, correlation=corPagel(1, kmtreedat$phy), weights=vf, method='ML')
gls13<- gls(Comp.2 ~ Comp.1 + np                                           , data=kmtreedat$data, correlation=corPagel(1, kmtreedat$phy), weights=vf, method='ML')
gls14<- gls(Comp.2 ~ Comp.1      + subf                                    , data=kmtreedat$data, correlation=corPagel(1, kmtreedat$phy), weights=vf, method='ML')
gls15<- gls(Comp.2 ~          np + subf                                    , data=kmtreedat$data, correlation=corPagel(1, kmtreedat$phy), weights=vf, method='ML')
gls16<- gls(Comp.2 ~ Comp.1                                                , data=kmtreedat$data, correlation=corPagel(1, kmtreedat$phy), weights=vf, method='ML')
gls17<- gls(Comp.2 ~          np                                           , data=kmtreedat$data, correlation=corPagel(1, kmtreedat$phy), weights=vf, method='ML')
gls18<- gls(Comp.2 ~               subf                                    , data=kmtreedat$data, correlation=corPagel(1, kmtreedat$phy), weights=vf, method='ML')
gls19<- gls(Comp.2 ~ 1                                                     , data=kmtreedat$data, correlation=corPagel(1, kmtreedat$phy), weights=vf, method='ML')

a<-sapply(list(gls1, gls2 ,gls3 ,gls4 ,gls5 ,gls6 ,gls7 ,gls8 ,gls9 ,gls10 ,gls11 ,gls12 ,gls13 ,gls14 ,gls15 ,gls16 ,gls17 ,gls18 ,gls19), function(x) AIC(x))
names(a)<- paste('gls',1:19,sep='')
sort(a, decreasing=F)
gls7      gls3      gls2      gls6      gls4      gls1      gls8      gls9     gls13      gls5     gls16     gls12     gls10     gls17     gls14     gls15 
-257.9115 -257.6655 -256.6656 -256.2829 -255.7958 -254.9979 -254.0335 -253.3195 -252.7938 -252.6640 -251.5338 -251.0484 -250.1295 -249.9888 -249.6610 -248.1910 
gls19     gls11     gls18 
-248.1335 -247.3184 -246.2370 

## Both top models contain np:Comp.1 interaction (gls7, 3). Nearest model without Comp.1:np is delta AIC~4.03
## Comp.1:subf interaction appears in gls3

### Next: work out akaike weights for these interactions.


## Look within each subfamily -- pattern holds within Eristalinae only
# Eristalinae

eris.phy1 <- drop.tip(kmtreedat$phy, tip=row.names(subset(kmtreedat$data, subf!='Eristalinae' | !np %in% c(4,5,6))))	
eris.dat1 <- within((a <- subset(kmtreedat$data, subf=='Eristalinae' & np %in% c(4,5,6))), binomial <- row.names(a))
eris.dat2 <- eris.dat1[,-grep('kcomp', names(eris.dat1))]
eris.dat3 <- subset(eris.dat1, !is.na(kcomp1)) 

kmtd.eris <- comparative.data(eris.phy1, eris.dat2, 'binomial')
kmtd.eris$data$genus <- factor(kmtd.eris$data$genus)
#		vf <- varFixed(~np)  
#		vf <- Initialize(vf, kmtd.eris$data)
gls1 <- gls(Comp.2 ~ SEX + Comp.1 * np, data=kmtd.eris$data, correlation=corPagel(1, kmtd.eris$phy), method='ML') # 2014-11-14
gls2 <- gls(Comp.2 ~ SEX + Comp.1 + np, data=kmtd.eris$data, correlation=corPagel(1, kmtd.eris$phy), method='ML')
anova(gls1, gls2)
#	     Model df       AIC       BIC   logLik   Test  L.Ratio p-value
#	gls1     1  7 -208.4034 -190.3804 111.2017                        
#	gls2     2  6 -203.9779 -188.5297 107.9890 1 vs 2 6.425424  0.0112	
with(kmtd.eris$data, xyplot(Comp.2 ~ Comp.1 | np))

summary(gls1) # for Supplementary table  # 2014-11-14
#	Generalized least squares fit by maximum likelihood
#	  Model: Comp.2 ~ SEX + Comp.1 * np 
#	  Data: kmtd.eris$data 
#	        AIC       BIC   logLik
#	  -208.4034 -190.3804 111.2017
#	
#	Correlation Structure: corPagel
#	 Formula: ~1 
#	 Parameter estimate(s):
#	   lambda 
#	0.9772449 
#	
#	Coefficients:
#	                  Value  Std.Error    t-value p-value
#	(Intercept)  0.12424088 0.17370016  0.7152606  0.4763
#	SEXM        -0.00574728 0.01477984 -0.3888596  0.6983
#	Comp.1      -0.26971101 0.11846768 -2.2766633  0.0251
#	np          -0.02783829 0.01291256 -2.1559068  0.0337
#	Comp.1:np    0.05342320 0.02126751  2.5119628  0.0137
#	
#	 Correlation: 
#	          (Intr) SEXM   Comp.1 np    
#	SEXM       0.067                     
#	Comp.1     0.092  0.133              
#	np        -0.354 -0.272 -0.324       
#	Comp.1:np -0.096 -0.153 -0.977  0.332
#	
#	Standardized residuals:
#	         Min           Q1          Med           Q3          Max 
#	-2.608551841 -0.250104873  0.004212361  0.403368257  0.825759392 
#	
#	Residual standard error: 0.2662667 
#	Degrees of freedom: 97 total; 92 residual

#### Use Head width for size
gls1 <- gls(Comp.2 ~ SEX + HW * np, data=kmtd.eris$data, correlation=corPagel(1, kmtd.eris$phy), method='ML') # 2014-11-14
gls2 <- gls(Comp.2 ~ SEX + HW + np, data=kmtd.eris$data, correlation=corPagel(1, kmtd.eris$phy), method='ML')
anova(gls1, gls2)


### kcomp1 as covariate
kmtd.eris <- comparative.data(eris.phy1, eris.dat3, 'binomial')
gls1 <- gls(Comp.2 ~ SEX + Comp.1 * kcomp1, data=kmtd.eris$data, correlation=corPagel(1, kmtd.eris$phy), method='ML') # 2014-11-14
gls2 <- gls(Comp.2 ~ SEX + Comp.1 + kcomp1, data=kmtd.eris$data, correlation=corPagel(1, kmtd.eris$phy), method='ML')
anova(gls1, gls2)
#	     Model df       AIC       BIC   logLik   Test  L.Ratio p-value
#	gls1     1  7 -77.73150 -65.08486 45.86575                        
#	gls2     2  6 -77.56894 -66.72897 44.78447 1 vs 2 2.162559  0.1414
summary(gls1)
#	Generalized least squares fit by maximum likelihood
#	  Model: Comp.2 ~ SEX + Comp.1 * kcomp1 
#	  Data: kmtd.eris$data 
#	       AIC       BIC   logLik
#	  -77.7315 -65.08486 45.86575
#	
#	Correlation Structure: corPagel
#	 Formula: ~1 
#	 Parameter estimate(s):
#	   lambda 
#	0.9473927 
#	
#	Coefficients:
#	                   Value  Std.Error    t-value p-value
#	(Intercept)   -0.1327217 0.09653525 -1.3748520  0.1768
#	SEXM          -0.0109626 0.02856297 -0.3838031  0.7032
#	Comp.1        -0.0010183 0.06321995 -0.0161079  0.9872
#	kcomp1         0.1787796 0.11314972  1.5800267  0.1220
#	Comp.1:kcomp1 -0.3718752 0.26495468 -1.4035426  0.1682
#	
#	 Correlation: 
#	              (Intr) SEXM   Comp.1 kcomp1
#	SEXM          -0.116                     
#	Comp.1        -0.064  0.111              
#	kcomp1         0.118  0.080 -0.529       
#	Comp.1:kcomp1 -0.148  0.297  0.591 -0.627
#	
#	Standardized residuals:
#	       Min         Q1        Med         Q3        Max 
#	-2.9044190  0.1645163  0.4843703  0.5918134  1.1914823 
#	
#	Residual standard error: 0.2011437 
#	Degrees of freedom: 45 total; 40 residual

### kcomp1 as covariate; HW as size measure
gls1 <- gls(Comp.2 ~ SEX + HW * kcomp1, data=kmtd.eris$data, correlation=corPagel(1, kmtd.eris$phy), method='ML') # 2014-11-14
gls2 <- gls(Comp.2 ~ SEX + HW + kcomp1, data=kmtd.eris$data, correlation=corPagel(1, kmtd.eris$phy), method='ML')
#		anova(gls1, gls2)
#	     Model df       AIC       BIC   logLik   Test  L.Ratio p-value
#	gls1     1  7 -77.64132 -64.99468 45.82066                        
#	gls2     2  6 -76.90707 -66.06709 44.45353 1 vs 2 2.734252  0.0982

summary(gls1)
#	Generalized least squares fit by maximum likelihood
#	  Model: Comp.2 ~ SEX + HW * kcomp1 
#	  Data: kmtd.eris$data 
#	        AIC       BIC   logLik
#	  -77.64132 -64.99468 45.82066
#	
#	Correlation Structure: corPagel
#	 Formula: ~1 
#	 Parameter estimate(s):
#	   lambda 
#	0.9406146 
#	
#	Coefficients:
#	                  Value  Std.Error    t-value p-value
#	(Intercept) -0.16121013 0.11433700 -1.4099559  0.1663
#	SEXM        -0.01172702 0.02903932 -0.4038326  0.6885
#	HW           0.01145733 0.02616898  0.4378209  0.6639
#	kcomp1      -0.29940770 0.26104185 -1.1469720  0.2582
#	HW:kcomp1    0.17154480 0.10831929  1.5836958  0.1211
#	
#	 Correlation: 
#	          (Intr) SEXM   HW     kcomp1
#	SEXM      -0.042                     
#	HW        -0.572 -0.096              
#	kcomp1     0.177  0.406 -0.467       
#	HW:kcomp1 -0.229 -0.307  0.577 -0.940
#	
#	Standardized residuals:
#	       Min         Q1        Med         Q3        Max 
#	-3.0126115  0.1276404  0.4984657  0.6322869  1.1691135 
#	
#	Residual standard error: 0.1958305 
#	Degrees of freedom: 45 total; 40 residual



## Using only subset with both PPCA data - 104 spp - is karyotype PPCA1 associated with allometry?  YES.
with(kmtreedat$data, xyplot(Comp.2 ~ Comp.1 | cut(kcomp1, breaks=fivenum(kcomp1))))
a<-within(kmtreedat$data, binomial<-row.names(kmtreedat$data))

kmtd.kcomp <- comparative.data(phy=drop.tip(kmtreedat$phy, tip=row.names(b<-subset(a, is.na(kcomp1)))),data=subset(a, !is.na(kcomp1), select=c('binomial','Comp.1','Comp.2','kcomp1','kcomp2','kcomp3','HW')),'binomial')
pgls1 <- pgls(Comp.2 ~ Comp.1 * kcomp1, data=kmtd.kcomp, lambda='ML') 	 ## 2014-09-15 
pgls2 <- pgls(Comp.2 ~ Comp.1 + kcomp1, data=kmtd.kcomp, lambda=pgls1$param[2])
anova(pgls1, pgls2)
#	Analysis of Variance Table
#	pgls: lambda = 0.96, delta = 1.00, kappa = 1.00
#	
#	Model 1: Comp.2 ~ Comp.1 * kcomp1
#	Model 2: Comp.2 ~ Comp.1 + kcomp1
#	  Res.Df    RSS Df Sum of Sq      F  Pr(>F)  
#	1    100 9.3183                              
#	2    101 9.8849 -1   -0.5666 6.0804 0.01537 *
#	---
#	Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#	
plot(resid(pgls1)~fitted(pgls1))  ## lovely
lrt(pgls1,pgls2)
$LR
[1] 6.138852

$DF
[1] 1

$P
[1] 0.01322431

summary(pgls1) # For Supplementary table
#	Call:
#	pgls(formula = Comp.2 ~ Comp.1 * kcomp1, data = kmtd.kcomp, lambda = "ML")
#	
#	Residuals:
#	     Min       1Q   Median       3Q      Max 
#	-0.62612 -0.22457 -0.01456  0.21668  0.75197 
#	
#	Branch length transformations:
#	
#	kappa  [Fix]  : 1.000
#	lambda [ ML]  : 0.956
#	   lower bound : 0.000, p = 1.8874e-15
#	   upper bound : 1.000, p = 1.2746e-06
#	   95.0% CI   : (0.863, 0.988)
#	delta  [Fix]  : 1.000
#	
#	Coefficients:
#	               Estimate Std. Error t value Pr(>|t|)   
#	(Intercept)   -0.045186   0.162206 -0.2786 0.781148   
#	Comp.1        -0.027843   0.041391 -0.6727 0.502696   
#	kcomp1         0.201408   0.073790  2.7295 0.007497 **
#	Comp.1:kcomp1 -0.392362   0.159118 -2.4659 0.015371 * 
#	---
#	Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#	
#	Residual standard error: 0.3053 on 100 degrees of freedom
#	Multiple R-squared: 0.09965,	Adjusted R-squared: 0.07264 
#	F-statistic: 3.689 on 3 and 100 DF,  p-value: 0.01443 

## Use Head width for size
pgls1 <- pgls(Comp.2 ~ HW * kcomp1, data=kmtd.kcomp, lambda='ML') 	 ## 2014-09-15 
pgls2 <- pgls(Comp.2 ~ HW + kcomp1, data=kmtd.kcomp, lambda=pgls1$param[2])
anova(pgls1, pgls2)
#	Analysis of Variance Table
#	pgls: lambda = 0.96, delta = 1.00, kappa = 1.00
#	
#	Model 1: Comp.2 ~ HW * kcomp1
#	Model 2: Comp.2 ~ HW + kcomp1
#	  Res.Df     RSS Df Sum of Sq      F   Pr(>F)   
#	1    100  9.6357                                
#	2    101 10.3356 -1  -0.69992 7.2638 0.008255 **

lrt(pgls1,pgls2)
#	$LR
#	[1] 7.292594
#	
#	$DF
#	[1] 1
#	
#	$P
#	[1] 0.006923942

summary(pgls1) # For Supplementary table
#	Call:
#	pgls(formula = Comp.2 ~ HW * kcomp1, data = kmtd.kcomp, lambda = "ML")
#	
#	Residuals:
#	    Min      1Q  Median      3Q     Max 
#	-0.7969 -0.1483 -0.0028  0.1897  0.8379 
#	
#	Branch length transformations:
#	
#	kappa  [Fix]  : 1.000
#	lambda [ ML]  : 0.961
#	   lower bound : 0.000, p = 1.1102e-15
#	   upper bound : 1.000, p = 1.8396e-06
#	   95.0% CI   : (0.876, 0.989)
#	delta  [Fix]  : 1.000
#	
#	Coefficients:
#	             Estimate Std. Error t value Pr(>|t|)   
#	(Intercept) -0.111514   0.171659 -0.6496 0.517425   
#	HW           0.023768   0.016678  1.4251 0.157242   
#	kcomp1      -0.286751   0.179635 -1.5963 0.113577   
#	HW:kcomp1    0.165030   0.061232  2.6951 0.008255 **
#	---
#	Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#	
#	Residual standard error: 0.3104 on 100 degrees of freedom
#	Multiple R-squared: 0.126,	Adjusted R-squared: 0.09978 
#	F-statistic: 4.806 on 3 and 100 DF,  p-value: 0.003613 



## Try again with heterogeneity
vf <- varFixed(~kcomp1) 
vf <- Initialize(vf, kmtd.kcomp$data)
gls1 <- gls(Comp.2 ~ Comp.1 * kcomp1, data=kmtd.kcomp$data, correlation=corPagel(1, kmtd.kcomp$phy),  method='ML')
gls2 <- gls(Comp.2 ~ Comp.1 + kcomp1, data=kmtd.kcomp$data, correlation=corPagel(1, kmtd.kcomp$phy),  method='ML')
gls1w <- gls(Comp.2 ~ Comp.1 * kcomp1, data=kmtd.kcomp$data, correlation=corPagel(1, kmtd.kcomp$phy),  weights=vf, method='ML')
gls2w <- gls(Comp.2 ~ Comp.1 + kcomp1, data=kmtd.kcomp$data, correlation=corPagel(1, kmtd.kcomp$phy),  weights=vf, method='ML')

anova(gls1,gls1w)
Model df        AIC       BIC   logLik
gls1      1  6 -170.31789 -154.4515 91.15895
gls1w     2  6  -56.42014  -40.5538 34.21007
### unweighted model is better by miles

### is karyotype PPCA2 associated with allometry?  NO.
pgls1 <- pgls(Comp.2 ~ Comp.1 * kcomp2, data=kmtd.kcomp, lambda='ML')   	 ## 2014-09-04 
pgls2 <- pgls(Comp.2 ~ Comp.1 + kcomp2, data=kmtd.kcomp, lambda=pgls1$param[2])
anova(pgls1,pgls2)
Analysis of Variance Table
pgls: lambda = 0.92, delta = 1.00, kappa = 1.00

Model 1: Comp.2 ~ Comp.1 * kcomp2
Model 2: Comp.2 ~ Comp.1 + kcomp2
Res.Df    RSS Df  Sum of Sq     F Pr(>F)
1    100 7.1813                           
2    101 7.1853 -1 -0.0040183 0.056 0.8135
> 
  
  gls1 <- gls(Comp.2 ~ Comp.1 * kcomp2, data=kmtd.kcomp$data, correlation=corPagel(1, kmtd.kcomp$phy),  method='ML')
gls2 <- gls(Comp.2 ~ Comp.1 + kcomp2, data=kmtd.kcomp$data, correlation=corPagel(1, kmtd.kcomp$phy),  method='ML')
anova(gls1, gls2)
Model df       AIC       BIC   logLik   Test    L.Ratio p-value
gls1     1  6 -162.7437 -146.8773 87.37183                          
gls2     2  5 -164.6855 -151.4635 87.34274 1 vs 2 0.05816797  0.8094
lrt(gls1,gls2)

plot(resid(pgls1)~fitted(pgls1))


