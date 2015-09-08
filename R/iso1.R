xixireg.isotonic <- function(x,y,k=3,nk=5){
    ds=data.frame(x=x,y=y)
    plot(x=ds$y,y=ds$x,col=1,type='l')
    isoreg(x=ds$y,y=ds$x) -> re.isoreg
    uug = diff(re.isoreg$yf)
    uug.inx = which(abs(uug)>1e-13)
    dds = data.frame(x = ds$x[-1][uug.inx], dfx = abs(diff(ds$y)[uug.inx]/(uug[uug.inx]) ))
    Kmeans(dds$dfx,3,iter.max=50,method = "manhattan") -> isoreg.knots
    Mode <- function(x) {
        ux <- unique(x)
        ux[which.max(tabulate(match(x, ux)))]
    }
    ggxfindk <- function(x,k=10){
        start = k
        n = length(x)
        re = NULL
        if (n < (100+2*k)) stop("x too short")
        for ( i in (k+1):(n-100-k)){
            m.pr = Mode(x[(i-k):i])
            m.after = Mode(x[i:(i+k)])
            if (m.pr != m.after) re = c(re,i)
        }
        return(re);
    }
    return(dds$x[ggxfindk(isoreg.knots$cluster,nk)]);
}
