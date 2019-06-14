## Computes the lower bound to the posterior expected Variation of Information

VI.lb=function(cls,psm){

	if(any(psm !=t(psm)) | any(psm >1) | any(psm < 0) | sum(diag(psm)) != nrow(psm) ){
	stop("psm must be a symmetric matrix with entries between 0 and 1 and 1's on the diagonals")}
    
	if(is.vector(cls)) cls <- t(cls)

	n=dim(psm)[1]

	VI.lb.compute=function(c){
		f=0
		for(i in 1:n){
			ind=(c==c[i])
			f=f+(log2(sum(ind))+log2(sum(psm[i,]))-2*log2(sum(ind*psm[i,])))/n
		}
		return(f)
	}
	output=apply(cls,1,VI.lb.compute)
	return(output)
}
