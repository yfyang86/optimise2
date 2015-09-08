#' Fix missing Column names.
#' 
#' @param x A matrix.
#' @return x with name pattern XIMP[i]
#' @examples
#' fix_missing_names(matrix(1:4,2))
fix_missing_names <-function(x){
    colnames(x) -> u
    if (!(is.null(u))){
        sum(u=="") -> n
        colnames(x)[u==""] = paste("XIMP",1:n,sep='')
    }else{
        n = dim(x)[2]
        colnames(x) = paste("XIMP",1:n,sep='')
    }
    x
}

cudaglmnetwrap <- function(x,y,lam){
    library(cudaglmnet)
	U = list();
	for (lami in lam) U[[as.character(lami)]] = cudaglmnet(x,y,lami)
	re = list(lambda = lam, beta = sapply(U,"[[","beta"))
	return(re);
}

#' Select Lambda using Dr Zhou's approach
#'
#' @param Yvec Y vector(n)
#' @param Xmat Design Matrix(n by p)
#' @param lambdas Proposed lambda, default value is seq(from=.001,to=1,length.out=50)
#' @param nPsep Number of independent simulation data.
#' @param perc Stopping percentage
#' @param iidSampletype Random number generator. Especially, users could define a custom function "r_usr_f".
#' @return result Result from glmnet
#' @return selection T/F of lambdas. F means
#' @examples
#' set.seed(65535)
#' Xmat = matrix(rnorm(100*80),ncol=80)
#' beta0 = rnorm(80,sd=2)
#' beta0[sample(1:80,70)] = 0.
#' epsilon = rnorm(100)
#' Yvec = Xmat%*%beta0 + epsilon
#' lambdas = seq(from=.001,to=1,length.out=50) # a vec
#' lambda_Select(Yvec,Xmat)
lambda_Select <- function(
                        Yvec,Xmat,lambdas=seq(from=.001,to=1,length.out=50),nPsep=20,perc = .2,
                        iidSampletype = "rnorm",cuda=FALSE,...){
    u = dim(Xmat)
    lambdas = sort(lambdas)
    # debug if (!("sizei"%in%ls())) sizei = 10L
	sizei=nPsep
   
    sampleN = u[1]
    samplep = u[2]
    if (!(iidSampletype%in%c('rnorm','rf','rt','rchisq','rbin'))){
    	iidSampletype = "rnorm";
    }else{
     if(iidSampletype == 'r_usr_f'){
     	iidSampletype=tryCatch({
     		 get(iidSampletype)(1)
     		 get(iidSampletype)(1:2)
     		 get(iidSampletype)(10)
     	},
     	error = function(err){
     	return("rnorm");
     	},
     	warning =  function(war){
     	return("rnorm");
     	}
     	)
     }
    }
    iidSampletype.f = get(iidSampletype)
    Zmat = matrix(iidSampletype.f(nPsep*sampleN),ncol=nPsep)
    uuz = paste("LPS",1:nPsep,sep='')
    colnames(Zmat) = uuz

    Xmat = fix_missing_names(cbind(Xmat,Zmat))
    if (!cuda){
            re = glmnet( Xmat,
            Yvec,
            lambda=lambdas,
            ...
            )}else{
			re = cudaglmnetwrap( 
			Xmat,
            Yvec,
            lambdas
            )
			
			}

    xlabel = colnames(Xmat)
    nperc = ceiling(sizei*perc)

    matchcol <- function(x){
        ifelse(sum(xlabel[abs(x)>(1e-14)] %in% uuz)> nperc,
               FALSE,
               TRUE
             )
     }
     re[[length(re)+1]] = perc
	   names(re)[length(re)] = "perc"
     re[[length(re)+1]] =samplep
	   names(re)[length(re)] = "samplep"
	   re[[length(re)+1]] = Yvec-Xmat%*%re$beta
     names(re)[length(re)] = "residual"
     re[[length(re)+1]] = apply(re$beta,2,function(x) sum(xlabel[abs(x)>(1e-14)] %in% uuz ))
	   names(re)[length(re)] = "L0"
     U = list(result=re,selection=apply(re$beta,2,matchcol));
     "lambda_select" -> class(U)
     return(U);
 }
 
 print.lambda_select <- function(x){
     cat('------------------\n')
     cat('Result from GLMNET\n')
     cat('------------------\n')
        print(x$result);
     cat('------------------\n')
     cat('Lambda Information\n Cut-off = ',
         x$result$perc *100,"%\n")
     cat('------------------\n')
        u = data.frame(Lambda = sort(x$result$lambda),
                       Select = x$selection,
                       NFalse = x$result[[length(x$result)]])
                       
        rownames(u) = paste("Lam",1:(dim(u)[1]),sep='')
        print(u,digits=5)
 }


summary.lambda_select <- function(x){
        u = data.frame(Lambda = sort(x$result$lambda),
                       Select = x$selection,
                       NFalse = x$result[[length(x$result)]])
                       
        rownames(u) = paste("Lam",1:(dim(u)[1]),sep='')
return(u); 
}


predict.lambda_select <- function(data,model){
  if (class(model)!="lambda_select") {warning("model MUST be a lambda_select class.");return(NULL);}
  data = as.matrix(data)
  data%*%model$result$beta[1:(model$result$samplep),]
}

