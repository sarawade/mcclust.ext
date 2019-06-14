### Calculate Credible ball

credibleball=function(c.star,cls.draw,c.dist=c("VI","Binder"),alpha=0.05){

	n=length(c.star)
	c.dist <- match.arg(c.dist, choices=c.dist)

	#Distance functions
	dist.binder=function(c1,c2){
		f=0
		for(i in 1:n){
			f=f+sum(abs((c1==c1[i])-(c2==c2[i])))
		}
		f=f/(n^2)
		return(f)
	}

	dist.vi=function(c1,c2){
		f=0
		for(i in 1:n){
			ind1=(c1==c1[i])
			ind2=(c2==c2[i])
			f=f+(log2(sum(ind1))+log2(sum(ind2))-2*log2(sum(ind1*ind2)))/n
		}
		return(f)
	}
	
	#Compute distance between optimal and samples
	M=nrow(cls.draw)
	d=rep(0,M)
	if(c.dist=="Binder") d=apply(cls.draw,1,dist.binder,c2=c.star)
	if(c.dist=="VI") d=apply(cls.draw,1,dist.vi,c2=c.star)
	sd=sort(d,decreasing=F,index.return=T)
	ind.star=ceiling((1-alpha)*M)
	
	cb=cls.draw[sd$ix[1:ind.star],]
	cb.dist=sd$x[1:ind.star]

	# Extremes of credible ball
	c.horiz=matrix(cb[(cb.dist==cb.dist[ind.star]),],ncol=n)
	k.cb=apply(cb,1,max)
	min.ind=which(k.cb==min(k.cb))
	c.uppervert=matrix(cb[min.ind[cb.dist[min.ind]==cb.dist[min.ind[length(min.ind)]]],],ncol=n)
	max.ind=which(k.cb==max(k.cb))
	c.lowervert=matrix(cb[max.ind[cb.dist[max.ind]==cb.dist[max.ind[length(max.ind)]]],],ncol=n)

	output=list(c.star=c.star,c.horiz=c.horiz,c.uppervert=c.uppervert,c.lowervert=c.lowervert,dist.horiz=cb.dist[ind.star],dist.uppervert=cb.dist[min.ind[length(min.ind)]],dist.lowervert=cb.dist[max.ind[length(max.ind)]])
	class(output)="credibleball"
	return(output)
}
