### Finds optimal partition minimizes the posterior expected Binder's loss

minbinder.ext=function(psm, cls.draw=NULL, method=c("avg","comp","draws","laugreen","greedy","all"), max.k=NULL, include.lg=FALSE, include.greedy=FALSE, start.cl.lg=NULL,start.cl.greedy=NULL,tol=0.001, maxiter=NULL,l=NULL, suppress.comment=TRUE){

	if(any(psm !=t(psm)) | any(psm >1) | any(psm < 0) | sum(diag(psm)) != nrow(psm) ){
	stop("psm must be a symmetric matrix with entries between 0 and 1 and 1's on the diagonals")}   	

	method <- match.arg(method, choices=method)
	if(method %in% c("draws","all") & is.null(cls.draw)) stop("cls.draw must be provided if method=''draws''")
	
	if(method != "greedy"){
		res=minbinder(psm, cls.draw,method=method,max.k=max.k,include.lg=include.lg,start.cl=start.cl.lg, tol=tol)
		n=nrow(psm)
		res$value=res$value/(n^2)*2
		if(method!="all"| (method=="all"&!include.greedy)) {
			if(method=="all") res=c(res,list(method="all"))
			class(res)="c.estimate"
			return(res)
		}
	} 

	if(method=="all" & is.null(start.cl.greedy)) {
		ind=which.min(res$value)
		start.cl.greedy=res$cl[ind,]
	}
      res.greedy <- greedy(psm,loss="Binder",start.cl=start.cl.greedy, maxiter=maxiter,L=l,suppress.comment=suppress.comment)
      if(method=="greedy") {
		res.greedy=c(res.greedy,list(method="greedy"))
		class(res.greedy)="c.estimate"
		return(res.greedy)  
	}

    	res$value <- c(res$value, res.greedy$value)
    	res$cl <- rbind(res$cl,res.greedy$cl)
    	res$cl[1,] <-res$cl[which.min(res$value),]
    	res$value[1] <- min(res$value)
	res=c(res,list(method="all"))
	rownames(res$cl)[nrow(res$cl)] <- names(res$value)[length(res$value)] <- "greedy"    
    	if(include.greedy) res=c(res,list(iter.greedy=res.greedy$iter.greedy))
    class(res)="c.estimate"
    return(res)     
}
