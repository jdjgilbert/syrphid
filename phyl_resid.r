phyl_resid<-function(C,x,Y) 
{ 
# find out how many y variables and taxa we have 
m<-ncol(Y); n<-nrow(Y); 
# compute inverse of C 
invC<-solve(C); 
# create a matrix for residuals 
r<-matrix(,n,m); 
# prepare X matrix 
X<-matrix(,n,2); X[,1]<-1; X[,2]<-x; 
# now loop over those variables, each time calculating the 
# regression & residuals 
for(i in 1:m) 
{ 
# estimate beta 
beta<- 
solve(t(X)%*%invC%*%X)%*%(t(X)%*%invC%*%Y[,i]); 
# compute residuals 
r[,i]<-Y[,i]-X%*%beta; 
} 
phyl_resid<-r; 
} 


