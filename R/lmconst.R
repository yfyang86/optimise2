lmconst<-function(y,x,x0,y0){
x=cbind(1,x)
n=dim(x)
p=x[2]
n=n[1]
beta=rep(0,)
A=matrix(c(1,x0),ncol=1)
b0=y0
Dmat=t(x)%*%x
dvec=t(x)%*%y
coef=solve.QP(Dmat,dvec,A,bvec=y0,meq=1)$solution
names(coef)[1]="Intercept"
coef
}

