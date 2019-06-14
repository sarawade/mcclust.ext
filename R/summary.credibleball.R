summary.credibleball=function(object,...){
	cb=object
	k.u=apply(cb$c.uppervert,1,max)
	n.u=nrow(cb$c.uppervert)
	k.l=apply(cb$c.lowervert,1,max)
	n.l=nrow(cb$c.lowervert)
	k.h=apply(cb$c.horiz,1,max)
	n.h=nrow(cb$c.horiz)

	t.u=list(n.u)
	for(i in 1:n.u){
		t=table(cb$c.star,cb$c.uppervert[i,],dnn=c("Estimate","Upper vertical bound"))
		t.u[[i]]=addmargins(t)
	}
	t.l=list(n.l)
	for(i in 1:n.l){
		t=table(cb$c.star,cb$c.lowervert[i,],dnn=c("Estimate","Lower vertical bound"))
		t.l[[i]]=addmargins(t)
	}
	t.h=list(n.h)
	for(i in 1:n.h){
		t=table(cb$c.star,cb$c.horiz[i,],dnn=c("Estimate","Horizontal bound"))
		t.h[[i]]=addmargins(t)
	}
	
	output=list(n.u=n.u,k.u=k.u,t.u=t.u,n.l=n.l,k.l=k.l,t.l=t.l,n.h=n.h,k.h=k.h,t.h=t.h,dist.uppervert=cb$dist.uppervert,dist.lowervert=cb$dist.lowervert,dist.horiz=cb$dist.horiz)
	class(output)="summary.credibleball"
	return(output)
}

print.summary.credibleball=function(x,...){
	cat("The credible ball characterizes the uncertainty in the clustering esitmate")
	cat(".\n")
	cat("It can be summarized with:")
	cat("\n")
 	cat("1. upper vertical bound: partitions in the ball with the fewest clusters that are most distant,")
	cat("\n")
      cat("2. lower vertical bound: partitions in the ball with the most clusters that are most distant,")
	cat("\n")
      cat("3. horizontal bound: partitions in the ball with the greatest distance.")
	cat("\n")	
	d.u=round(x$dist.uppervert,2)
	d.l=round(x$dist.lowervert,2)
	d.h=round(x$dist.horiz,2)
	if(x$n.u==1){
		if(x$k.u!=1) {
			cat(paste("The upper vertical bound has",x$k.u,"clusters with a distance of", d.u))
			cat(".\n")
		}
		else {
			cat(paste("The upper vertical bound has 1 cluster with a distance of", d.u))
			cat(".\n")
		}
	}
	else{
		cat(paste("The", x$n.u, "upper vertical bounds have", x$k.u[1], "clusters with a distance of", d.u))
		cat(".\n")
	}

	if(x$n.l==1){
		cat(paste("The lower vertical bound has", x$k.l, "clusters with a distance of", d.l))
		cat(".\n")
	}
	else{
		cat(paste("The", x$n.l, "lower vertical bounds have", x$k.l[1], "clusters with a distance of", d.l))
		cat(".\n")
	}

	if(x$n.h==1){
		cat(paste("The horizontal bound has", x$k.h, "clusters with a distance of", d.h))
		cat(".\n")
	}
	else{
		cat(paste("The", x$n.h, "horizontal bounds have "))
		for(i in 1:(x$n.h-1)){
			cat(x$k.h[i])
			if(i<(x$n.h-1)){ cat(", ")}
		}
		cat(paste(" and",x$k.h[x$n.h], "clusters with a distance of", d.h))
		cat(".\n")
	}

	cat("\nA cross tabulation with the point estimate of the partition and the bounds.\n")
	cat("\n")
	for (i in 1:x$n.u){
		print(x$t.u[[i]])
	}
	cat("\n")
	for (i in 1:x$n.l){
		print(x$t.l[[i]])
	}
	cat("\n")
	for (i in 1:x$n.h){
		print(x$t.h[[i]])
	}

}
