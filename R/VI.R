## Computes the posterior expected Variation of Information

VI=function(cls,cls.draw){

	if(is.vector(cls)) cls <- t(cls)

	n=dim(cls.draw)[2]
	M=dim(cls.draw)[1]

	VI.compute=function(c){
		f=0
		for(i in 1:n){
			ind=(c==c[i])
			f=f+log2(sum(ind))
			for(m in 1:M){
				indm=(cls.draw[m,]==cls.draw[m,i])
				f=f+(log2(sum(indm))-2*log2(sum(ind*indm)))/M
			}
		}
		f=f/n
		return(f)
	}
	output=apply(cls,1,VI.compute)
	return(output)
}
