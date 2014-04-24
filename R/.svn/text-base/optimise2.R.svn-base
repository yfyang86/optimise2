optimise2 <-
function(f, interval, ..., lower = min(interval), upper = max(interval), 
    maximum = FALSE,tol = .Machine$double.eps^0.25,initials=runif(1,lower,upper),trace=F,maxit=100){

optimise2R<-function (f, interval, ..., lower = min(interval), upper = max(interval), 
    maximum = FALSE,tol = .Machine$double.eps^0.25,initials=runif(1,lower,upper),trace=T){
    trace0=trace+0;
    if (maximum) {
val <- .External2("optimise2", function(arg) -f(arg, ...), 
    lower, upper, tol,initials,trace0,package="optimise2")
list(maximum = val, objective = f(val, ...))
    }
    else {
val <- .External2("optimise2", function(arg) f(arg, ...), 
    lower, upper, tol,initials,trace0,package="optimise2")
list(minimum = val, objective = f(val, ...))
    }
}

x=matrix(0,ncol=2,nrow=maxit);
for (i in 1:maxit) {
re=optimise2R(f, c(1.5, 5),trace=F)
x[i,]=c(re[[1]],re[[2]]);
}
if (maximum) {id=which.max(x[,2])[1]}else{id =which.min(x[,2])[1]}
val=x[id,1];
 
if (maximum) return(list(maximum = val, objective = f(val, ...)));
return(list(minimum = val, objective = f(val, ...)));
}
