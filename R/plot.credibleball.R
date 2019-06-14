plot.credibleball=function(x,data=NULL,dx=NULL,xgrid=NULL,dxgrid=NULL,...){
	if(is.null(data)){
		stop("data must be supplied")}
	if(!is.data.frame(data)){
		stop("data must be a data.frame")}
	p=ncol(data)
	n=nrow(data)
	k.u=apply(x$c.uppervert,1,max)
	n.u=nrow(x$c.uppervert)
	k.l=apply(x$c.lowervert,1,max)
	n.l=nrow(x$c.lowervert)
	k.h=apply(x$c.horiz,1,max)
	n.h=nrow(x$c.horiz)
 	par(ask = TRUE)
	if(p==1){
		x1=data[,1]
		if(is.null(dxgrid)){
			if(is.null(xgrid)){
				dout=density(x1)
			}
			else{
				dout=density(data,n=length(xgrid),from=xgrid[1],to=xgrid[length(xgrid)])
			}
			xgrid=dout$x
			dxgrid=dout$y
		}
		#Compute density estimate at x1
		x2=dx
		if(is.null(x2)){
			aout=approx(xgrid,dxgrid,xout=x1)
			x2=aout$y}

		#Upper
		cl=rainbow(k.u)
		par(mfrow=c(1,n.u),mar=c(2.5,2.5,1.5,.5)+.1, mgp=c(1.5, .5, 0))
		for(i in 1:n.u){	
			plot(xgrid,dxgrid,type="l",xlab=names(data),ylab="density",main="Credible Ball: Upper Vertical Bound",...)
			for(j in 1:k.u){
				points(x1[x$c.uppervert[i,]==j], x2[x$c.uppervert[i,]==j],col=cl[j],...)
		}}
		#Lower
		cl=rainbow(k.l)
		par(mfrow=c(1,n.l),mar=c(2.5,2.5,1.5,.5)+.1, mgp=c(1.5, .5, 0))
		for(i in 1:n.l){	
			plot(xgrid,dxgrid,type="l",xlab=names(data),ylab="density",main="Credible Ball: Lower Vertical Bound",...)
			for(j in 1:k.l){
				points(x1[x$c.lowervert[i,]==j], x2[x$c.lowervert[i,]==j],col=cl[j],...)
		}}
		#Horiziontal
		par(mfrow=c(1,n.h),mar=c(2.5,2.5,1.5,.5)+.1, mgp=c(1.5, .5, 0))
		for(i in 1:n.h){
			cl=rainbow(k.h[i])	
			plot(xgrid,dxgrid,type="l",xlab=names(data),ylab="density",main="Credible Ball: Horizontal Bound",...)
			for(j in 1:k.h[i]){
				points(x1[x$c.horiz[i,]==j], x2[x$c.horiz[i,]==j],col=cl[j],...)
		}}
	}
	if(p==2){
		x1=data[,1]
		x2=data[,2]
		#Upper
		cl=rainbow(k.u)
		par(mfrow=c(1,n.u),mar=c(2.5,2.5,1.5,.5)+.1, mgp=c(1.5, .5, 0))
		for(i in 1:n.u){
			plot(x1,x2,xlab=names(data)[1],ylab=names(data)[2],main="Credible Ball: Upper Vertical Bound",...)
			for(j in 1:k.u){
				points(x1[x$c.uppervert[i,]==j], x2[x$c.uppervert[i,]==j],col=cl[j],...)
		}}
		#Lower
		cl=rainbow(k.l)
		par(mfrow=c(1,n.l),mar=c(2.5,2.5,1.5,.5)+.1, mgp=c(1.5, .5, 0))
		for(i in 1:n.l){	
			plot(x1,x2,xlab=names(data)[1],ylab=names(data)[2],main="Credible Ball: Lower Vertical Bound",...)
			for(j in 1:k.l){
				points(x1[x$c.lowervert[i,]==j], x2[x$c.lowervert[i,]==j],col=cl[j],...)
		}}
		#Horiziontal
		par(mfrow=c(1,n.h),mar=c(2.5,2.5,1.5,.5)+.1, mgp=c(1.5, .5, 0))
		for(i in 1:n.h){
			cl=rainbow(k.h[i])	
			plot(x1,x2,xlab=names(data)[1],ylab=names(data)[2],main="Credible Ball: Horizontal Bound",...)
			for(j in 1:k.h[i]){
				points(x1[x$c.horiz[i,]==j], x2[x$c.horiz[i,]==j],col=cl[j],...)
		}}
	}
	if(p>2){
		x.pca=princomp(data,scores=T)
		x1=x.pca$scores[,1]
		x2=x.pca$scores[,2]
		#Upper
		cl=rainbow(k.u)
		par(mfrow=c(1,n.u),mar=c(2.5,2.5,1.5,.5)+.1, mgp=c(1.5, .5, 0))
		for(i in 1:n.u){
			plot(x1,x2,xlab="PC 1",ylab="PC 2",main="Credible Ball: Upper Vertical Bound",...)
			for(j in 1:k.u){
				points(x1[x$c.uppervert[i,]==j], x2[x$c.uppervert[i,]==j],col=cl[j],...)
		}}
		#Lower
		cl=rainbow(k.l)
		par(mfrow=c(1,n.l),mar=c(2.5,2.5,1.5,.5)+.1, mgp=c(1.5, .5, 0))
		for(i in 1:n.l){	
			plot(x1,x2,xlab="PC 1",ylab="PC 2",main="Credible Ball: Lower Vertical Bound",...)
			for(j in 1:k.l){
				points(x1[x$c.lowervert[i,]==j], x2[x$c.lowervert[i,]==j],col=cl[j],...)
		}}
		#Horiziontal
		par(mfrow=c(1,n.h),mar=c(2.5,2.5,1.5,.5)+.1, mgp=c(1.5, .5, 0))
		for(i in 1:n.h){
			cl=rainbow(k.h[i])	
			plot(x1,x2,xlab="PC 1",ylab="PC 2",main="Credible Ball: Horizontal Bound",...)
			for(j in 1:k.h[i]){
				points(x1[x$c.horiz[i,]==j], x2[x$c.horiz[i,]==j],col=cl[j],...)
		}}
	}	
}
