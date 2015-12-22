

#### Analysis of morphology WRT different taxonomic levels

mtreedat$data$genus <- factor(sub('_.*','',row.names(mtreedat$data)))
mtreedat$data$tri <- syrph_tax$tri[match(mtreedat$data$genus, syrph_tax$genus)]
mtreedat$data$subtri <- syrph_tax$subtri[match(mtreedat$data$genus, syrph_tax$genus)]
mtreedat$data$subf <- syrph_tax$subf[match(mtreedat$data$tri, syrph_tax$tri)]

pgls.subf <- pgls(Comp.2 ~ Comp.1 * subf + SEX - 1, data=mtreedat, lambda='ML')   	 ## 2014-09-04 
pgls.subf1 <- pgls(Comp.2 ~ Comp.1 + subf + SEX - 1, data=mtreedat, lambda=pgls.subf$param[2])
anova(pgls.subf, pgls.subf1)
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

lrt(pgls.subf, pgls.subf1)
> lrt(pgls.subf, pgls.subf1)
$LR
[1] 9.834297

$DF
[1] 1

$P
[1] 0.001712877
