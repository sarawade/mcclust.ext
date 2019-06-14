plot.c.estimate=function(x,data=NULL,dx=NULL,xgrid=NULL,dxgrid=NULL,...){
	if(is.null(data)){
		stop("data must be supplied")}
	if(!is.data.frame(data)){
		stop("data must be a data.frame")}
	p=ncol(data)
	n=nrow(data)
	if(!is.matrix(x$cl)) x$cl=matrix(x$cl,nrow=1)
	k=apply(x$cl,1,max)
	n.c=nrow(x$cl)
 	par(ask = TRUE)
	method.names=x$method
	if(method.names=="all"){method.names=names(x$value)}
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
		for(i in 1:n.c){
			cl=rainbow(k[i])
			plot(xgrid,dxgrid,type="l",xlab=names(data),ylab="density",main=paste("Method:",method.names[i]),...)
			for(j in 1:k[i]){
				points(x1[x$cl[i,]==j], x2[x$cl[i,]==j],col=cl[j],...)
		}}

	}
	if(p==2){
		x1=data[,1]
		x2=data[,2]
		for(i in 1:n.c){
			cl=rainbow(k[i])
			plot(x1,x2,xlab=names(data)[1],ylab=names(data)[2],main=paste("Method:",method.names[i]),...)
			for(j in 1:k[i]){
				points(x1[x$cl[i,]==j], x2[x$cl[i,]==j],col=cl[j],...)
		}}
	}
	if(p>2){
		x.pca=princomp(data,scores=T)
		x1=x.pca$scores[,1]
		x2=x.pca$scores[,2]
		for(i in 1:n.c){
			cl=rainbow(k[i])
			plot(x1,x2,xlab="PC 1",ylab="PC 2",main=paste("Method:",method.names[i]),...)
			for(j in 1:k[i]){
				points(x1[x$cl[i,]==j], x2[x$cl[i,]==j],col=cl[j],...)
		}}
	}	
}
