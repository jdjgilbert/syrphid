#The function below performs distance-based phylogenetic least- squares analysis. The method may
#be used to assess regression, ANOVA, and other linear mod- els in a phylogenetic context.
#The approach is particularly useful for high-dimensional data where standard (parametric)
#phylogenetic generalized least squares analyses, or the anal- ysis of phylogenetically
#independent contrasts, cannot be performed. The function below obtains the sums of squares
#(SS) for all factors in the linear model, and statistically evaluates them via permutation.

D.pgls<-function(f1,phy,iter = 999){ 
library(ape)
data = NULL
form.in<-formula(f1)
Terms<-terms(form.in,keep.order = TRUE)
Y<-as.matrix(eval(form.in[[2]],parent.frame()))
N<-length(phy$tip.label)
p<-ncol(Y)
if(is.null(rownames(Y))){ stop("No species names with Y-data.")}
if(length(match(rownames(Y), phy$tip.label))!= N)
stop("Data matrix missing some taxa present on the tree.")
if(length(match(phy$tip.label,rownames(Y)))!= N)
stop("Tree missing some taxa in the data matrix.")
C<-vcv.phylo(phy)
C<-C[rownames(Y),rownames(Y)]
eigC <- eigen(C)
D.mat<-solve(eigC$vectors%*%
diag(sqrt(eigC$values))%*% t(eigC$vectors))
Y.new<-D.mat%*% (Y)
ones.new<-D.mat%*%(array(1,N))
pred.1<- predict(lm(Y.new~ones.new-1))
dat<-model.frame(form.in,data)
df<-df.tmp<-SS.tmp<-SS.obs<-F<-array()
for (i in 1:ncol(attr(Terms, "factors"))){
mod.mat<-model.matrix(Terms[1:i],data = dat)
x.new<-D.mat%*%mod.mat
pred.y<-predict(lm(Y.new~x.new-1))
G<-(pred.y-pred.1)%*%t(pred.y-pred.1)
SS.tmp[i]<-sum(diag(G))
ifelse(i == 1, SS.obs[i]<-SS.tmp[i], SS.obs[i]<- (SS.tmp[i]-SS.tmp[i-1]))
df.tmp[i]<-ifelse(ncol(mod.mat) == 1,1,(ncol(mod.mat)- 1))
ifelse(i == 1, df[i]<-df.tmp[i], df[i]<-(df.tmp[i]-df.tmp[i- 1]))
}
MS<-SS.obs/df
mod.mat<-model.matrix(Terms)
x.new<-D.mat%*%mod.mat
y.res<-residuals(lm(Y.new~x.new-1))
SS.res<-sum(diag(y.res%*%t(y.res)))
df.res<-nrow(Y)-1-sum(df)
MS.res<-SS.res/df.res
Rsq<-SS.obs/(sum(SS.obs)+SS.res)
F<-MS/MS.res
F.r<-P.val<-array(1,dim = length(SS.obs))
for(i in 1:iter){
SS.tmp<-SS.r<-array()
Y.r<-as.matrix(Y[sample(nrow(Y)),])
row.names(Y.r)<-row.names(Y)
 Y.r.new<-D.mat%*% (Y.r)
pred.1.r<- predict(lm(Y.r.new~ones.new-1))
for (ii in 1:ncol(attr(Terms, "factors"))){
mod.mat<-model.matrix(Terms[1:ii])
x.new<-D.mat%*%mod.mat
pred.y.r<-predict(lm(Y.r.new~x.new-1))
G.r<-(pred.y.r-pred.1.r)%*%t(pred.y.r-pred.1.r)
SS.tmp[ii]<-sum(diag(G.r))
ifelse(ii == 1, SS.r[ii]<-SS.tmp[ii], SS.r[ii]<-(SS.tmp[ii]- SS.tmp[ii-1]))
}
MS.r<-SS.r/df
mod.mat<-model.matrix(Terms)
x.new<-D.mat%*%mod.mat
y.res.r<-residuals(lm(Y.r.new~x.new-1))
SS.r.res<-sum(diag(y.res.r%*%t(y.res.r)))
MS.r.res<-SS.r.res/df.res
F.r<-MS.r/MS.r.res
P.val<-ifelse(F.r>= F, P.val+1,P.val)
}
P.val<-P.val/(iter+1)
anova.tab<-cbind(df,SS.obs,MS,F,P.val,Rsq)
anova.tab<-rbind(anova.tab,c(df.res,SS.res,MS.res,NA,NA,NA))
rownames(anova.tab)<-c(colnames(attr(Terms, "factors")), "Residual")
return(anova.tab)
}



dpgls.data <- kmtreedat$data






