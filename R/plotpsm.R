plotpsm=function(psm,method="complete",...){

	if(any(psm !=t(psm)) | any(psm >1) | any(psm < 0) | sum(diag(psm)) != nrow(psm) ){
	stop("psm must be a symmetric matrix with entries between 0 and 1 and 1's on the diagonals")}
    
	#sort by heirarchical clustering
	hc=hclust(as.dist(1-psm), method = method, members = NULL)
	psm_hc=psm
	n=nrow(psm)
	psm_hc[1:n,]=psm_hc[hc$order,]
	psm_hc[,1:n]=psm_hc[,hc$order]

	image(1:n,1:n,1-psm_hc,col=heat.colors(20),...)
}

