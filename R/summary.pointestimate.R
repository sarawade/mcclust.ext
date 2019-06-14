summary.c.estimate=function(object,...){
	x=object
	if(!is.matrix(x$cl)) x$cl=matrix(x$cl,nrow=1)
	k=apply(x$cl,1,max)
	n.c=nrow(x$cl)

	t=list(n.c)
	onevec=rep("size",ncol(x$cl))
	for(i in 1:n.c){
		t[[i]]=table(onevec,x$cl[i,],dnn=c("","cluster"))
	}

	output=list(method=x$method,k=k,n.c=n.c,t=t,value=x$value)
	class(output)="summary.c.estimate"
	return (output)
}

print.summary.c.estimate=function(x,...){

	if(x$method!="all"){
		cat("The partition estimate found with the",x$method, "method has a posterior expected loss of\n", round(x$value,2),"and contains",x$k,"clusters of sizes:\n")
		print(x$t[[1]])
	}
	if(x$method=="all"){
		cat("The best partition estimate has a posterior expected loss of\n", round(x$value[1],2)," and contains",x$k[1],"clusters of sizes:\n")
		print(x$t[[1]])
		cat("\n")
		method.names=names(x$value)
		for(i in 2:x$n.c){
			cat("The partition estimate found with the",method.names[i], "method has a posterior expected loss of\n", round(x$value[i],2)," and contains",x$k[i],"clusters of sizes:\n")
			print(x$t[[i]])
			cat("\n")
		}
	}
}
