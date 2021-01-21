lslr2<-function(x, dist="weibull", npar=2, reg_method="XonY")  {
	## a convergence limit is fixed here for 3rd parameter  convergence
	## no longer an argument for the R function, but still an argument to C++ functions
	limit<-1e-5

	if(is.vector(x))  {
		stop("use MRR functions for casual fitting, or pre-process with getPPP")
	}else{
		if(names(x)[1]=="time"&&names(x)[2]=="ppp")  {
		## will handle the output from getPPP
		}else{
			if(length(x$ppp)<3)  {
				stop("insufficient failure points")
			}else{
				stop("input format not recognized")
			}
		}
	}

## It turns out that this code is general to all fitting methods:
	if(tolower(dist) %in% c("weibull","weibull2p","weibull3p")){
		fit_dist<-"weibull"
	}else{
		if(tolower(dist) %in% c("lnorm", "lognormal","lognormal2p", "lognormal3p")){
			fit_dist<-"lnorm"
		}else{
			if(!dist=="gumbel") {
		## Note: lslr contains experimental support for "gumbel"
			stop(paste0("dist argument ", dist, "is not recognized for distribution fitting"))
			}
		}
	}

##	npar<-2 ## introducing 3p in dist argument will override any npar (or its default)
	if(tolower(dist) %in% c("weibull3p", "lognormal3p")){
		npar<-3
	}

	##casenum<-0
	reg_order=0
	if(reg_method=="YonX") reg_order=1
	##if(npar==3) casenum=casenum+2
	dist_num=0
	if(fit_dist=="lnorm")dist_num=1
	if(dist=="gumbel") dist_num=3

	parlist<-list(fail=x$time, ppp=x$ppp, reg_order=reg_order, dist_num=dist_num, npar=npar, limit=limit)
	result<-.Call("LSLR2", parlist , package="ParaBoot")
	return(result)
	
	## Note this is not a complete rework of original lslr
	## continued code would assign names to the output vector elements
	## and also call for and provide the AbPval, which really should be optional
	## so that it can be elimintated for call from mlefit for the start vector.
}
