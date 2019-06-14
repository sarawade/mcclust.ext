### Compute optimal partition that minimizes the posterior expected loss,
### by performing greedy search: at every iteration, consider the L-closest 
### ancestors and the L-closest descendents

greedy=function(psm,cls.draw=NULL,loss=NULL,start.cl=NULL,maxiter=NULL,L=NULL,suppress.comment=TRUE){

	if(any(psm !=t(psm)) | any(psm >1) | any(psm < 0) | sum(diag(psm)) != nrow(psm) ){
	stop("psm must be a symmetric matrix with entries between 0 and 1 and 1's on the diagonals")}

	n=nrow(psm)
	if(is.null(loss)) loss="VI.lb"
	if(is.null(start.cl)) start.cl=1:n
	if(is.null(maxiter)) maxiter=2*n
	if(is.null(L)) L=2*n

	if(loss=="VI" & is.null(cls.draw)) stop("cls.draw must be provided if loss=''VI''")

	EVI_lb_local=function(c){
		f=0
		for(i in 1:n){
			ind=(c==c[i])
			f=f+(log2(sum(ind))+log2(sum(psm[i,]))-2*log2(sum(ind*psm[i,])))/n
		}
		return(f)
	}
	EVI_local=function(c){
		M=nrow(cls.draw)
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
	EBL_local=function(c){
		f=0
		for(i in 1:n){
			f=f+sum(abs((c[i]==c)-psm[i,]))
		}
		f=f/(n^2)
		return(f)
	}  

	#Extra functions

	c_combine=function(c,i,j){
		c[c==i|c==j]=min(i,j)
		c[c>max(i,j)]=c[c>max(i,j)]-1
		return(c)
	}

	dist_merge_ij=function(ni,nj){
		d=0
		if(loss=="VI.lb"||loss=="VI"){d=((ni+nj)/n)*log2((ni+nj)/n)-(ni/n)*log2(ni/n)-(nj/n)*log2(nj/n)}
		if(loss=="Binder"){d=((ni+nj)^2-(ni)^2-(nj)^2)/(n^2)}
		return(d)
	}
	dist_split_i=function(x,ni){
		d=0
		if(loss=="VI.lb"||loss=="VI"){d=(ni/n)*log2(ni/n)-(x/n)*log2(x/n)-((ni-x)/n)*log2((ni-x)/n)}
		if(loss=="Binder"){d=((ni)^2-(x)^2-(ni-x)^2)/(n^2)}
		return(d)
	}

	#Function which given a configuration, finds the L closests configurations and
	# selects the one with the smallest EBL
	local_explore=function(c_star,val_star){
		k=max(c_star)
		nj=rep(0,k)
		for(j in 1:k){
			nj[j]=sum(c_star==j)
		}
		snj_ind=list()
		unj=unique(nj)
		unj=sort(unj)
		U=length(unj)
		lnj=rep(0,U)
		for(i in 1:U){
			snj_ind[[i]]=which(nj==unj[i])
			lnj[i]=length(snj_ind[[i]])
		}
		c_opt=c_star
		val_opt=val_star
		#Merge two clusters
		#Compute distance of merge any two clusters
		if(k>1){
			m_ind=1:U
			if(lnj[1]==1){m_ind=m_ind[-1]}
			d_1=apply(matrix(unj[m_ind],length(m_ind),1),1,dist_merge_ij,nj=unj[1])
			d_mat=rbind(d_1,rep(1,length(m_ind)),m_ind)
			if((U-(lnj[U]==1))>1){
				for(i in 2:(U-(lnj[U]==1))){
					m_ind=i:U
					if(lnj[i]==1){m_ind=m_ind[-1]}
					d_i=apply(matrix(unj[m_ind],length(m_ind),1),1,dist_merge_ij,nj=unj[i])
					d_mat=cbind(d_mat,rbind(d_i,rep(i,length(m_ind)),m_ind))
				}
			}
			sd=sort(d_mat[1,],index.return=T)
			d_mat=matrix(d_mat[,sd$ix],nrow=3)
			colind=1
			l=0
			ub=min(L,choose(k,2))
			while(l<ub){
				i=d_mat[2,colind]
				h=d_mat[3,colind]
				reps=0
				if(i!=h){reps=lnj[i]*lnj[h]}
				if(i==h){reps=choose(lnj[i],2)}
				nj_indi=snj_ind[[i]]
				if(i==h){nj_indi=nj_indi[-lnj[i]]}
				for(i_ind in 1:length(nj_indi)){
					n1_ind=nj_indi[i_ind]
					if(i!=h){nj_indh=snj_ind[[h]]}
					if(i==h){nj_indh=snj_ind[[i]][(i_ind+1):lnj[i]]}
					for(h_ind in 1:length(nj_indh)){
						n2_ind=nj_indh[h_ind]
						#proposed partition
						c_p=c_combine(c_star,n1_ind,n2_ind)
						#compute loss
						if(loss=="VI.lb"){val_p=EVI_lb_local(c_p)}
						if(loss=="VI"){val_p=EVI_local(c_p)}
						if(loss=="Binder"){val_p=EBL_local(c_p)}
						if(val_p<val_opt){
							c_opt=c_p
							val_opt=val_p
						}
					}
				}
				#Update l and colind
				colind=colind+1
				l=l+reps
			}
		}
		#Spliting two clusters
		#Compute distance of splitting any clusters
		if(k<n){
			sind=1+(unj[1]==1)
			m_ind=1:floor(unj[sind]/2)
			d_1=apply(matrix(m_ind,length(m_ind),1),1,dist_split_i,ni=unj[sind])
			d_mat=rbind(d_1,rep(sind,length(m_ind)),m_ind)
			numsp=apply(matrix(m_ind,length(m_ind),1),1,choose,n=unj[sind])
			if((unj[sind]%%2)==0){numsp[length(numsp)]=numsp[length(numsp)]/2}
			numsp=sum(numsp)*lnj[sind]
			if(sind<U){
				for(i in (sind+1):U){
					m_ind=1:floor(unj[i]/2)
					d_i=apply(matrix(m_ind,length(m_ind),1),1,dist_split_i,ni=unj[i])
					d_mat=cbind(d_mat,rbind(d_i,rep(i,length(m_ind)),m_ind))
					numsp=c(numsp,apply(matrix(m_ind,length(m_ind),1),1,choose,n=unj[i]))
					if((unj[i]%%2)==0){numsp[length(numsp)]=numsp[length(numsp)]/2}
					numsp=numsp[1]+sum(numsp[-1])*lnj[i]
				}
			}
			sd=sort(d_mat[1,],index.return=T)
			d_mat=matrix(d_mat[,sd$ix],nrow=3)
			colind=1
			l=0
			ub=min(L,numsp)
			while(l<ub){
				i=d_mat[2,colind]
				nj_new=d_mat[3,colind]
				reps=choose(unj[i],nj_new)
				if(nj_new==(unj[i]/2)){reps=reps/2}
				for(j in 1:lnj[i]){
					ind_set=c(1:nj_new)
					for(h in 1:reps){
						c_p=c_star
						c_p[c_star==snj_ind[[i]][j]][ind_set]=k+1
						#Compute expected loss
						val_p=0
						if(loss=="VI.lb"){val_p=EVI_lb_local(c_p)}
						if(loss=="VI"){val_p=EVI_local(c_p)}
						if(loss=="Binder"){val_p=EBL_local(c_p)}
						if(val_p<val_opt){
							c_opt=c_p
							val_opt=val_p
						}
						if(h<reps){
							#Update set
							ind_set[nj_new]=ind_set[nj_new]+1
							if(ind_set[nj_new]>unj[i]){
								updateind=which(ind_set>=c((unj[i]-nj_new+1):unj[i]))[1]
								ind_set[c((updateind-1):nj_new)]=ind_set[(updateind-1)]+c(1:(nj_new-updateind+2))
							}
						}
					}
				}
				colind=colind+1
				l=l+reps*lnj[i]	

			}
		}
		return(list(c_star=c_opt,val_star=val_opt))
	}

	#Start at last
	c_star=start.cl
	val_star=0
	if(loss=="VI.lb") val_star=EVI_lb_local(c_star)
	if(loss=="VI") val_star=EVI_local(c_star)
	if(loss=="Binder") val_star=EBL_local(c_star)

	it=1
	stop_ind=F
	while((it<maxiter)&(!stop_ind)){
		opt=local_explore(c_star,val_star)
		if(opt$val_star==val_star){
			stop_ind=T
		}
		else{
			val_star=opt$val_star
			c_star=opt$c_star
			it=it+1
			if(!suppress.comment){cat(paste("Iteration=",it," k=",max(c_star), " Loss=",round(val_star,4),"\n"))}
		}
	}

output=list(cl=c_star,value=val_star,iter.greedy=it)
return(output)
}
