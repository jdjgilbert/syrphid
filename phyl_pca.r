phyl_pca<-function(C,X,mode) { 
# find out how many columns and taxa we have 
m<-ncol(X); n<-nrow(X); 
# compute inverse of C 
invC<-solve(C); 
# compute vector of ancestral states 
one<-matrix(1,n,1); 
a<-t(t(one)%*%invC%*%X)*sum(sum(invC))^-1; 
# compute evolutionary VCV matrix 
V<- t(X-one%*%t(a))%*%invC%*%(X-one%*%t(a))*(n-1)^-1; 
# if correlation matrix 
if(mode == 'corr'){ 
# standardize X 
X = X/(one%*%t(sqrt(diag(V)))); 
# change V to correlation matrix 
V = V/(sqrt(diag(V))%*%t(sqrt(diag(V)))); 
# recalculate a 
a<-t(t(one)%*%invC%*%X)*sum(sum(invC))^-1; 
} 
# eigenanalyze 
es = eigen(V); 
result<-NULL; result$Eval<-diag(es$values); 
result$Evec <-es$vectors; 
# compute scores in the species space 
result$S<-(X-one%*%t(a))%*%result$Evec; 
# compute cross covariance matrix 
# and loadings 
Ccv<-t(X-one%*%t(a))%*%invC%*%result$S/(n-1); 
result$L<-matrix(,m,m); 
for(i in 1:m) 
for(j in 1:m) 
result$L[i,j]<-Ccv[i,j]/sqrt(V[i,i]*result$Eval[j,j]); 
phyl_pca<-result; 
# done 
}
