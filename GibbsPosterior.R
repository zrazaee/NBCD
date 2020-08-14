GibbsPosterior<- function(N,alpha,beta,n,z,nSamp, burnT){
  K <-  prod(N)
  Tmax <- burnT + nSamp
  
  # generating the directed lattice
  A <- matrix(0,K,K)
  for (i in 1:(K-1)){
    if  (i%%N[2]!=0) {A[i,i+1]= 1}
    if  (i<= (N[2]*(N[1]-1))){
      A[i,i+N[2]]=1
    }
  }

  lattice_sample = latticeGibbsSamplerC(A, alpha+z, beta+n-z, nSamp, burnT)
  lattice_sample = lattice_sample[complete.cases(lattice_sample), ]
  median = apply(lattice_sample,2,median)  
  sd = apply(lattice_sample,2,sd)
  mean = apply(lattice_sample,2,mean)
  return(list(median=median, mean = mean, sd=sd,lattice_sample=lattice_sample))
}

