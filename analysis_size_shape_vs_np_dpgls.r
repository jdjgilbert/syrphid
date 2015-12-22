
### D-PGLS analysis of size-shape relationship WRT NP (or KPPCA1)

### Same question using DPGLS
row.names(kpdat) <- kpdat$binomial
temp <- kpdat[!is.na(kmtreedat$data$np),grep("WL", names(kpdat)):grep("PS", names(kpdat))]
Y <- as.matrix(temp)		### Y is a matrix of all morphological data with np data as well
#row.names(Y) <- row.names(kmtreedat$data)
Y <- Y[match(row.names(kmtreedat$data), row.names(Y)),]

### MPPCA1 as size measure

procD.pgls(Y ~ kmtreedat$data$SEX + kmtreedat$data$Comp.1 + kmtreedat$data$np + kmtreedat$data$Comp.1:kmtreedat$data$np, phy=kmtreedat$phy, iter=1000)
#                                         df      SS      MS     Rsq        F   P.val
#kmtreedat$data$SEX                        1   40514   40514 0.00328   2.1115 0.45255
#kmtreedat$data$Comp.1                     1 8870530 8870530 0.71882 462.3115 0.00100
#kmtreedat$data$np                         1   32173   32173 0.00261   1.6768 0.27772
#kmtreedat$data$Comp.1:kmtreedat$data$np   1  384848  384848 0.03119  20.0574 0.00100
#Residuals                               157 3012413   19187      

### HW as size measure
temp <- kpdat[!is.na(kmtreedat$data$np),grep("WL", names(kpdat)):grep("PS", names(kpdat))]
X <- grep("HW", names(temp))
Y <- as.matrix(temp[,-X])		### Y is a matrix of all morphological data with np data as well
Y <- Y[match(row.names(kmtreedat$data), row.names(Y)),]

procD.pgls(Y ~ kmtreedat$data$SEX + kmtreedat$data$HW + kmtreedat$data$np + kmtreedat$data$HW:kmtreedat$data$np, phy=kmtreedat$phy, iter=1000)
#	                                     df      SS      MS     Rsq        F   P.val
#	kmtreedat$data$SEX                    1   40374   40374 0.00327   2.1653 0.49251
#	kmtreedat$data$HW                     1 8864568 8864568 0.71873 475.4165 0.00100
#	kmtreedat$data$np                     1   10229   10229 0.00083   0.5486 0.55644
#	kmtreedat$data$HW:kmtreedat$data$np   1  491155  491155 0.03982  26.3412 0.00100
#	Residuals                           157 2927406   18646     




#### Use DPGLS incorporating subfamily


procD.pgls(Y ~ kmtreedat$data$Comp.1 + kmtreedat$data$np + kmtreedat$data$subf + kmtreedat$data$Comp.1:kmtreedat$data$np + kmtreedat$data$Comp.1:kmtreedat$data$subf, phy=kmtreedat$phy, iter=1000)                                      df      SS      MS     Rsq        F   P.val
kmtreedat$data$Comp.1                       1 8854071 8854071 0.71748 499.9116 0.00100
kmtreedat$data$np                           1   23270   23270 0.00189   1.3138 0.31369
kmtreedat$data$subf                         1      29      29 0.00000   0.0016 0.55245
kmtreedat$data$Comp.1:kmtreedat$data$np     1  353860  353860 0.02867  19.9794 0.00200				##... but NP still sig effect
kmtreedat$data$Comp.1:kmtreedat$data$subf   1  346290  346290 0.02806  19.5519 0.00100  ## subfamilies have different allometries
Residuals                                 156 2762958   17711 

procD.pgls(Y ~ kmtreedat$data$Comp.1 * kmtreedat$data$np * kmtreedat$data$subf, phy=kmtreedat$phy, iter=1000)
df      SS      MS     Rsq        F   P.val
kmtreedat$data$Comp.1                                         1 8854071 8854071 0.71748 500.5043 0.00100
kmtreedat$data$np                                             1   23270   23270 0.00189   1.3154 0.29471
kmtreedat$data$Comp.1:kmtreedat$data$np                       1  353880  353880 0.02868  20.0042 0.00100
kmtreedat$data$subf                                           1       8       8 0.00000   0.0005 0.80420
kmtreedat$data$Comp.1:kmtreedat$data$subf                     1  346290  346290 0.02806  19.5751 0.00100
kmtreedat$data$np:kmtreedat$data$subf                         1   36424   36424 0.00295   2.0590 0.12687
kmtreedat$data$Comp.1:kmtreedat$data$np:kmtreedat$data$subf   1    2228    2228 0.00018   0.1259 0.84216 ## no evidence for 3-way interaction
Residuals                                                   154 2724306   17690 



procD.pgls(Y ~ kmtreedat$data$SEX + kmtreedat$data$Comp.1 + kmtreedat$data$np + kmtreedat$data$genus + kmtreedat$data$Comp.1:kmtreedat$data$np + kmtreedat$data$Comp.1:kmtreedat$data$genus, phy=kmtreedat$phy, iter=500)
#	procD.pgls(Y ~ kmtreedat$data$Comp.1 + kmtreedat$data$np + kmtreedat$data$genus + kmtreedat$data$Comp.1:kmtreedat$data$np + kmtreedat$data$Comp.1:kmtreedat$data$genus, phy=kmtreedat$phy, iter=500)
#	                                           df      SS      MS     Rsq        F    P.val
#	kmtreedat$data$Comp.1                       1 8854071 8854071 0.71748 761.0230 0.001996
#	kmtreedat$data$np                           1   23270   23270 0.00189   2.0001 0.053892
#	kmtreedat$data$genus                       66  499791    7573 0.04050   0.6509 0.001996
#	kmtreedat$data$Comp.1:kmtreedat$data$np     1  362896  362896 0.02941  31.1916 0.001996
#	kmtreedat$data$Comp.1:kmtreedat$data$genus 66 2297955   34817 0.18621   2.9926 0.001996
#	Residuals                                  26  302495   11634   

## drop np interaction
> procD.pgls(Y ~ kmtreedat$data$SEX + kmtreedat$data$Comp.1 + kmtreedat$data$np + kmtreedat$data$genus + kmtreedat$data$Comp.1:kmtreedat$data$genus, phy=kmtreedat$phy, iter=500)
df      SS      MS     Rsq        F    P.val
kmtreedat$data$SEX                          1   40514   40514 0.00328   4.0088 0.109780
kmtreedat$data$Comp.1                       1 8870530 8870530 0.71882 877.7267 0.001996
kmtreedat$data$np                           1   32173   32173 0.00261   3.1835 0.025948
kmtreedat$data$genus                       66  475306    7202 0.03852   0.7126 0.001996
kmtreedat$data$Comp.1:kmtreedat$data$genus 66 2659193   40291 0.21549   3.9867 0.001996
Residuals                                  26  262763   10106 

## drop genus interaction
> procD.pgls(Y ~ kmtreedat$data$SEX + kmtreedat$data$Comp.1 + kmtreedat$data$np + kmtreedat$data$genus + kmtreedat$data$Comp.1:kmtreedat$data$np, phy=kmtreedat$phy, iter=500)
df      SS      MS     Rsq        F   P.val
kmtreedat$data$SEX                       1   40514   40514 0.00328   1.4596 0.46307
kmtreedat$data$Comp.1                    1 8870530 8870530 0.71882 319.5801 0.00200
kmtreedat$data$np                        1   32173   32173 0.00261   1.1591 0.23154
kmtreedat$data$genus                    66  475306    7202 0.03852   0.2595 0.20359
kmtreedat$data$Comp.1:kmtreedat$data$np  1  396085  396085 0.03210  14.2698 0.00200
Residuals                               91 2525871   27757                         
> 
  
  
  ##### DPGLS, using kcomp1
  
  
  kmtd.kcomp <- subset(kmtreedat, !is.na(kcomp1)) 
Y <- as.matrix(subset(kpdat, row.names(kpdat)%in%row.names((subset(kmtd.kcomp$data,!is.na(kcomp1)))))[,grep("WL", names(kpdat)):grep("PS", names(kpdat))])		### Y is a matrix of all morphological data with np data as well
Y <- Y[match(row.names(kmtd.kcomp$data), row.names(Y)),]
procD.pgls(Y ~ kmtd.kcomp$data$SEX + kmtd.kcomp$data$Comp.1 * kmtd.kcomp$data$kcomp1, phy=kmtd.kcomp$phy, iter=1000)
#	                                              df      SS      MS     Rsq        F   P.val
#	kmtd.kcomp$data$SEX                            1   16739   16739 0.00702   3.2153 0.37163
#	kmtd.kcomp$data$Comp.1                         1 1729139 1729139 0.72551 332.1435 0.00100
#	kmtd.kcomp$data$kcomp1                         1   45017   45017 0.01889   8.6471 0.02198
#	kmtd.kcomp$data$Comp.1:kmtd.kcomp$data$kcomp1  1   77069   77069 0.03234  14.8040 0.00400
#	Residuals                                     99  515394    5206     

temp <- subset(kpdat, row.names(kpdat)%in%row.names(kmtd.kcomp$data))[,grep("WL", names(kpdat)):grep("PS", names(kpdat))]
X <- grep("HW", names(temp))
Y <- as.matrix(temp[,-X])		### Y is a matrix of all morphological data with np data as well
Y <- Y[match(row.names(kmtd.kcomp$data), row.names(Y)),]

procD.pgls(Y ~ kmtd.kcomp$data$SEX + kmtd.kcomp$data$HW * kmtd.kcomp$data$kcomp1, phy=kmtd.kcomp$phy, iter=1000)


#### Only Eristalinae, using DPGLS

#### MPPCA1 as size measure

Y <- as.matrix(kpdat[which(row.names(kpdat)%in%row.names(kmtd.eris$data)),grep("WL", names(kpdat)):grep("PS", names(kpdat))])		### Y is a matrix of all morphological data with np data as well
#Y <- Y[match(row.names(kmtd.eris$data), row.names(Y)),]

procD.pgls(Y ~ kmtd.eris$data$SEX + kmtd.eris$data$Comp.1 * kmtd.eris$data$np, phy=kmtd.eris$phy, iter=1000)
#	                                        df      SS      MS     Rsq        F   P.val
#	kmtd.eris$data$SEX                       1   46261   46261 0.00409   1.7836 0.50749
#	kmtd.eris$data$Comp.1                    1 8636106 8636106 0.76442 332.9734 0.00100
#	kmtd.eris$data$np                        1   62794   62794 0.00556   2.4211 0.26174
#	kmtd.eris$data$Comp.1:kmtd.eris$data$np  1  166217  166217 0.01471   6.4086 0.03896
#	Residuals                               92 2386142   25936 


### HW as size measure
temp <- subset(kpdat, !is.na(kmtreedat$data$np)&row.names(kpdat)%in%row.names(kmtd.eris$data))[,grep("WL", names(kpdat)):grep("PS", names(kpdat))]
X <- grep("HW", names(temp))
Y <- as.matrix(temp[,-X])		### Y is a matrix of all morphological data with np data as well
Y <- Y[match(row.names(kmtd.eris$data), row.names(Y)),]

procD.pgls(Y ~ kmtd.eris$data$SEX + kmtd.eris$data$HW + kmtd.eris$data$np + kmtd.eris$data$HW:kmtd.eris$data$np, phy=kmtd.eris$phy, iter=1000)
#	                                    df      SS      MS     Rsq        F    P.val
#	kmtd.eris$data$SEX                   1   15169   15169 0.00970   2.5495 0.245754
#	kmtd.eris$data$HW                    1 1163956 1163956 0.74401 195.6260 0.000999
#	kmtd.eris$data$np                    1   66067   66067 0.04223  11.1039 0.019980
#	kmtd.eris$data$HW:kmtd.eris$data$np  1   63404   63404 0.04053  10.6563 0.001998
#	Residuals                           43  255846    5950



### Kcomp1, MPPCA1 as size measure

Y <- as.matrix(subset(kpdat, row.names(kpdat)%in%row.names((subset(kmtd.eris$data,!is.na(kcomp1)))))[,grep("WL", names(kpdat)):grep("PS", names(kpdat))])		### Y is a matrix of all morphological data with np data as well

kmtd.eris.kc1 <- subset(kmtd.eris, subf=='Eristalinae')

procD.pgls(Y ~ kmtd.eris.kc1$data$SEX + kmtd.eris.kc1$data$Comp.1 + kmtd.eris.kc1$data$kcomp1 + kmtd.eris.kc1$data$Comp.1:kmtd.eris.kc1$data$kcomp1, phy=kmtd.eris.kc1$phy, iter=1000)

### Kcomp1, HW as size measure
Y <- as.matrix(subset(kpdat, row.names(kpdat)%in%row.names((subset(kmtd.eris$data,!is.na(kcomp1)))))[,grep("WL", names(kpdat)):grep("PS", names(kpdat))])		### Y is a matrix of all morphological data with np data as well
X <- grep("HW", names(temp))
Y <- as.matrix(temp[,-X])		### Y is a matrix of all morphological data with np data as well
Y <- Y[match(row.names(kmtd.eris$data), row.names(Y)),]

kmtd.eris.kc1 <- subset(kmtd.eris, subf=='Eristalinae')

procD.pgls(Y ~ kmtd.eris.kc1$data$SEX + kmtd.eris.kc1$data$HW + kmtd.eris.kc1$data$kcomp1 + kmtd.eris.kc1$data$HW:kmtd.eris.kc1$data$kcomp1, phy=kmtd.eris.kc1$phy, iter=1000)



kmtreedat$data$kcomp1 <- ktreedat$data$Comp.1[match(row.names(kmtreedat$data), row.names(ktreedat$data))]
kmtreedat$data$kcomp2 <- ktreedat$data$Comp.2[match(row.names(kmtreedat$data), row.names(ktreedat$data))]
kmtreedat$data$kcomp3 <- ktreedat$data$Comp.3[match(row.names(kmtreedat$data), row.names(ktreedat$data))]
kmtreedat$data$kcomp4 <- ktreedat$data$Comp.4[match(row.names(kmtreedat$data), row.names(ktreedat$data))]



