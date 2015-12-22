##########  Use D.pgls function from Adams 2014 to look at multivariate data ###########




#### dpgls.r is deprecated as it is now incorporated into package geomorph
#source('dpgls.r') ### Contains code for the D.pgls() function -- removes need for PPCA for shape

library(geomorph)

### We're looking at allometry, so retain the PPCA1 size axis
 Y <- as.matrix(mtreedat$data[,2:grep('PS', names(mtreedat$data)),])  ### Y is a matrix of all morphological data

### procD.pgls is the official version from library(geomorph)
procD.pgls(Y ~ mtreedat$data$Comp.1 * mtreedat$data$genus + mtreedat$data$SEX, phy=mtreedat$phy, iter=1000)
#	                                          df       SS       MS     Rsq         F    P.val
#	mtreedat$data$Comp.1                       1 18860010 18860010 0.73232 2426.1436 0.000999
#	mtreedat$data$genus                       82   775057     9452 0.03009    1.2159 0.001998
#	mtreedat$data$Comp.1:mtreedat$data$genus  82  4720842    57571 0.18331    7.4059 0.000999
#	mtreedat$data$SEX                          1   154226   154226 0.00599   19.8396 0.001998
#	Residuals                                160  1243785     7774  

### Use head width instead to represent body size
procD.pgls(Y[,-grep('HW', colnames(Y))] ~ mtreedat$data$HW * mtreedat$data$genus + mtreedat$data$SEX, phy=mtreedat$phy, iter=1000)
#	                                      df       SS       MS     Rsq         F    P.val
#	mtreedat$data$HW                       1 17963486 17963486 0.69801 2189.5723 0.000999
#	mtreedat$data$genus                   82   961922    11731 0.03738    1.4299 0.000999
#	mtreedat$data$HW:mtreedat$data$genus  82  5434244    66271 0.21116    8.0778 0.000999
#	mtreedat$data$SEX                      1    62795    62795 0.00244    7.6541 0.061938  ### Sex only marginally sig in this model
#	Residuals                            160  1312657     8204    

## subtribe
procD.pgls(Y ~ mtreedat$data$Comp.1 * mtreedat$data$subtri + mtreedat$data$SEX, phy=mtreedat$phy, iter=1000)
# df       SS       MS     Rsq         F    P.val
# mtreedat$data$Comp.1                        1 18860010 18860010 0.73232 1861.0340 0.000999
# mtreedat$data$subtri                       22   139205     6327 0.00541    0.6244 0.001998
# mtreedat$data$Comp.1:mtreedat$data$subtri  22  3722276   169194 0.14453   16.6955 0.000999
# mtreedat$data$SEX                           1   194866   194866 0.00757   19.2287 0.009990
# Residuals                                 280  2837564    10134  

procD.pgls(Y[,-grep('HW', colnames(Y))] ~ mtreedat$data$HW * mtreedat$data$subtri + mtreedat$data$SEX, phy=mtreedat$phy, iter=1000)
# df       SS       MS     Rsq         F    P.val
# mtreedat$data$HW                        1 17963486 17963486 0.69801 1696.8619 0.000999
# mtreedat$data$subtri                   22   113205     5146 0.00440    0.4861 0.004995
# mtreedat$data$HW:mtreedat$data$subtri  22  4589554   208616 0.17834   19.7062 0.000999
# mtreedat$data$SEX                       1   104695   104695 0.00407    9.8897 0.051948
# Residuals                             280  2964163    10586 

## tribe
procD.pgls(Y ~ mtreedat$data$Comp.1 * mtreedat$data$tri + mtreedat$data$SEX, phy=mtreedat$phy, iter=1000)
# df       SS       MS     Rsq         F    P.val
# mtreedat$data$Comp.1                     1 18860010 18860010 0.73232 1568.7143 0.000999
# mtreedat$data$tri                       13    15567     1197 0.00060    0.0996 0.137862
# mtreedat$data$Comp.1:mtreedat$data$tri  13  3136020   241232 0.12177   20.0649 0.000999
# mtreedat$data$SEX                        1   159591   159591 0.00620   13.2742 0.033966
# Residuals                              298  3582732    12023       

procD.pgls(Y[,-grep('HW', colnames(Y))] ~ mtreedat$data$HW * mtreedat$data$tri + mtreedat$data$SEX, phy=mtreedat$phy, iter=1000)
# df       SS       MS     Rsq         F    P.val
# mtreedat$data$HW                     1 17963486 17963486 0.69801 1274.7411 0.000999
# mtreedat$data$tri                   13    88254     6789 0.00343    0.4818 0.007992
# mtreedat$data$HW:mtreedat$data$tri  13  3432873   264067 0.13339   18.7390 0.000999
# mtreedat$data$SEX                    1    51114    51114 0.00199    3.6272 0.256743
# Residuals                          298  4199377    14092        


  
## Same thing for np
row.names(kpdat) <- kpdat$binomial
kpdat <- kpdat[match(row.names(kmtreedat$data), row.names(kpdat)),]
Y <- as.matrix(kpdat[,grep("WL", names(kpdat)):grep("PS", names(kpdat))])		### Y is a matrix of all morphological data with np data as well

procD.pgls(Y ~ kmtreedat$data$Comp.1 * kmtreedat$data$np + kmtreedat$data$SEX, phy=kmtreedat$phy, iter=1000)
#	                                         df      SS      MS     Rsq        F   P.val
#	kmtreedat$data$Comp.1                     1 8854071 8854071 0.71748 461.4537 0.00050
#	kmtreedat$data$np                         1   23270   23270 0.00189   1.2128 0.32784
#	kmtreedat$data$Comp.1:kmtreedat$data$np   1  353880  353880 0.02868  18.4434 0.00150
#	kmtreedat$data$SEX                        1   96844   96844 0.00785   5.0473 0.25987
#	Residuals                               157 3012413   19187      

kmtreedat$data$HW <- mtreedat$data$HW[match(row.names(kmtreedat$data), row.names(mtreedat$data))]
procD.pgls(Y[,-grep('HW', colnames(Y))] ~ kmtreedat$data$HW * kmtreedat$data$np + kmtreedat$data$SEX, phy=kmtreedat$phy, iter=1000)
#	                                     df      SS      MS     Rsq        F   P.val
#	kmtreedat$data$HW                     1 8891106 8891106 0.72088 476.8398 0.00100
#	kmtreedat$data$np                     1    8574    8574 0.00070   0.4599 0.56144
#	kmtreedat$data$HW:kmtreedat$data$np   1  468310  468310 0.03797  25.1160 0.00100
#	kmtreedat$data$SEX                    1   38335   38335 0.00311   2.0559 0.47253
#	Residuals                           157 2927406   18646    

Y <- as.matrix(kpdat[!is.na(kmtreedat$data$kcomp1),grep("WL", names(kpdat)):grep("PS", names(kpdat))])		### Y is a matrix of all morphological data with np data as well

### Look at KPPCA1
kmtreedat.kcomp1 <- subset(kmtreedat$data, !is.na(kcomp1))
kmtreedat.kcomp1$binomial <- row.names(kmtreedat.kcomp1)
kmtreedat.kcomp1$genus <- 
kmtreedat.kcomp1 <- comparative.data(drop.tip(kmtreedat$phy, tip=which(!kmtreedat$phy$tip.label %in% row.names(kmtreedat.kcomp1))), kmtreedat.kcomp1, binomial)

kmtreedat.kcomp1$data$genus <- do.call(rbind.data.frame, strsplit(row.names(kmtreedat.kcomp1$data), '_'))[,1]
#kmtreedat.kcomp1$data$genus <- factor(kmtreedat.kcomp1$data$genus)

procD.pgls(Y ~ kmtreedat.kcomp1$data$SEX + kmtreedat.kcomp1$data$Comp.1 + kmtreedat.kcomp1$data$kcomp1 + kmtreedat.kcomp1$data$genus + kmtreedat.kcomp1$data$Comp.1:kmtreedat.kcomp1$data$kcomp1, phy=kmtreedat.kcomp1$phy, iter=500)
#	                                                          df      SS      MS     Rsq        F    P.val
#	kmtreedat.kcomp1$data$SEX                                  1   16739   16739 0.00702   3.3205 0.301397
#	kmtreedat.kcomp1$data$Comp.1                               1 1729139 1729139 0.72551 343.0092 0.001996
#	kmtreedat.kcomp1$data$kcomp1                               1   45017   45017 0.01889   8.9300 0.007984
#	kmtreedat.kcomp1$data$genus                               47  284023    6043 0.11917   1.1988 0.003992
#	kmtreedat.kcomp1$data$Comp.1:kmtreedat.kcomp1$data$kcomp1  1   46303   46303 0.01943   9.1852 0.045908
#	Residuals                                                 52  262136    5041                          

  procD.pgls(Y ~ kmtreedat.kcomp1$data$SEX + kmtreedat.kcomp1$data$Comp.1 + kmtreedat.kcomp1$data$kcomp1 + kmtreedat.kcomp1$data$genus + kmtreedat.kcomp1$data$Comp.1:kmtreedat.kcomp1$data$genus, phy=kmtreedat.kcomp1$phy, iter=500)
#	                                                         df      SS      MS     Rsq        F    P.val
#	kmtreedat.kcomp1$data$SEX                                 1   16739   16739 0.00702   1.1108 0.187625
#	kmtreedat.kcomp1$data$Comp.1                              1 1729139 1729139 0.72551 114.7503 0.001996
#	kmtreedat.kcomp1$data$kcomp1                              1   45017   45017 0.01889   2.9875 0.001996
#	kmtreedat.kcomp1$data$genus                              47  284023    6043 0.11917   0.4010 0.001996
#	kmtreedat.kcomp1$data$Comp.1:kmtreedat.kcomp1$data$genus 47  218028    4639 0.09148   0.3078 0.065868
#	Residuals                                                 6   90412   15069                          
> 
  procD.pgls(Y ~ kmtreedat.kcomp1$data$SEX + kmtreedat.kcomp1$data$Comp.1 + kmtreedat.kcomp1$data$kcomp1 + kmtreedat.kcomp1$data$genus + kmtreedat.kcomp1$data$Comp.1:kmtreedat.kcomp1$data$kcomp1 + kmtreedat.kcomp1$data$Comp.1:kmtreedat.kcomp1$data$genus, phy=kmtreedat.kcomp1$phy, iter=500)
#                                                          df      SS      MS     Rsq       F    P.val
#kmtreedat.kcomp1$data$SEX                                  1   16739   16739 0.00702  0.9610 0.185629
#kmtreedat.kcomp1$data$Comp.1                               1 1729139 1729139 0.72551 99.2768 0.001996
#kmtreedat.kcomp1$data$kcomp1                               1   45017   45017 0.01889  2.5846 0.005988
#kmtreedat.kcomp1$data$genus                               47  284023    6043 0.11917  0.3470 0.001996
#kmtreedat.kcomp1$data$Comp.1:kmtreedat.kcomp1$data$kcomp1  1   46303   46303 0.01943  2.6585 0.015968
#kmtreedat.kcomp1$data$Comp.1:kmtreedat.kcomp1$data$genus  47  175050    3724 0.07345  0.2138 0.101796
#Residuals                                                  5   87087   17417                         
#> 

