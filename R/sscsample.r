sscsample<-function(size, n.samples, sample.type="simple", x = NULL, 
					strata = NULL, cluster = NULL, cl.size = NULL, ret=FALSE)
{
	# size - the sample size 
	# n.samples - the number of samples to draw
	# sample.type - the method of sampling can be one of :
	#      "cluster","stratified","simple"

	# DO NOT set the following parameters unless you know what you're doing!
	# x - a data vector
	# strata - the stratum each data point in x belongs to
	# cluster - the cluster each data point in x belongs to
	# cl.size - the number of clusters to sample 
	#	- must be less than the number of clusters

	# if ret is true then the samples and their summary statistics are
	# return in a list
	
	data(sscsample.data)

	if(is.null(x))
		x<-sscsample.data$value

	nx<-length(x)

	if(size>nx)
		stop("Sample size must be less than population size")

	if(is.null(strata))
		strata<-sscsample.data$stratum
	
	# just in case the strata numbers are not evenly spaced
	strata.names<-unique(strata)
	n.strata<-length(strata.names)

	if(nx!=length(strata))
		stop("The length of the strata and data vectors must be equal")

	if(is.null(cluster))
		cluster<-sscsample.data$cluster

	n.clusters<-length(unique(cluster))

	if(nx!=length(cluster))
		stop("The length of the cluster and data vectors must be equal")

	samples<-matrix(0,nrow=size,ncol=n.samples)

	for(r in 1:n.samples){
		if(sample.type=="simple"|sample.type==1){
			sample.idx<-sample(1:nx,size)

		} else if(sample.type=="stratified"|sample.type==2){
			
			stratified.data<-split(1:nx,strata.names)

			k<-size*sapply(stratified.data,length)/nx
			
			for(stratum in 1:n.strata){
				if(stratum==1)
					sample.idx<-sample(stratified.data[[stratum]],k[stratum])
				else
					sample.idx<-c(sample.idx,sample(stratified.data[[stratum]],k[stratum]))
			}	
		} else if(sample.type=="cluster"|sample.type==3){
			# just in case the cluster numbers are not evenly spaced			
			
			cluster.names<-unique(cluster) 

			if(cl.size<0)
				cl.size<-4

			if(cl.size>n.clusters)
				stop("The number of randomly sampled clusters must be less that the overall number of clusters")
		
			cluster.sample<-sample(cluster.names,cl.size)

			clustered.data<-split(1:nx,cluster.names)
			
			for(i in 1:cl.size)
			{
				cl<-cluster.sample[i]
				if(i==1)
					sample.idx<-clustered.data[[cl]]
				else
					sample.idx<-c(sample.idx,clustered.data[[cl]])
			}		
		} else
			stop(paste("Unknown sampling sample.type :",sample.type))
	
		samples[,r]<-sample.idx
	}

	means<-rep(0,n.samples)
	s.strata<-matrix(0,nrow=n.samples,ncol=n.strata)
	sample.out<-matrix(0,nrow=size,ncol=n.samples)

	cat("Sample\tMean   \tStratum 1\tStratum 2\tStratum 3\n")
	cat("------\t-------\t---------\t---------\t---------\n")

	for(r in 1:n.samples)
	{
		idx<-samples[,r]
		means[r]<-mean(x[idx])

		for(j in 1:n.strata)
			s.strata[r,j]<-sum(strata[idx]==strata.names[j])

		sample.out[,r]<-x[idx]

		cat(paste(r,"\t",round(means[r],4),"\t",s.strata[r,1],"\t\t",s.strata[r,2]
			,"\t\t",s.strata[r,3],"\n",sep=""))
	}
	
	if(ret)
		return(list(samples=samples,s.strata=s.strata,means=means))
	
}
