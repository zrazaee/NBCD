source('GibbsPosterior.R')
ncores= 15
library(LearnBayes)
d <- c(4,4)
K <- prod(d)
pmin <- 0.04
pmax <- 0.34
quantile1 <- list(p= 0.5^(1/K),x=pmax)
quantile2 <- list(p= 1-0.5^(1/K), x = pmin)
hyper_select <- beta.select(quantile1,quantile2)


A <- matrix(0,K,K)
for (i in 1:(K-1)){
  if  (i%%d[2]!=0) {A[i,i+1]= 1}
  if  (i<= (d[2]*(d[1]-1))){
    A[i,i+d[2]]= 1
  }
}
alpha <- beta <- rep(1,K)
Rcpp::sourceCpp('RcppArmaGibbs.cpp')



m <-  seq(hyper_select[1],sum(hyper_select),length.out= 15)
counter <- 0

gibbs_med <- gibbs_sd <- par <- c()

res = result = list()
hyperparameter_fx <- function(t, t2){
  x = seq(0.2,t/2, length.out= 10)
  y = seq(0.2,t2/2, length.out= 10)
  z = seq(0.2,0.4,length.out= 3)
  for (i in x){
    for (j in y){
      for (l in z){
        gibbs_sample = GibbsPosterior(d,c(t-i,rep(l,K-2),j),c(i,rep(min(t/2,t2/2)-l,K-2),t2-j),rep(0,K),rep(0,K), 5000, 1000)
        gibbs_med = rbind(gibbs_med,gibbs_sample$median)
        gibbs_sd = rbind(gibbs_sd,gibbs_sample$sd)
        par = rbind(par,c(t,t2,i,j,l))
        counter = counter+1
        cat(i)
        cat("\n")
      }
    }
    #}
  }
  res = list(gibbs_med=gibbs_med,gibbs_sd=gibbs_sd, par=par )
}

library(doParallel)
library(foreach)
comb <- function(...){
  mapply('rbind', ..., SIMPLIFY = FALSE)
}

registerDoParallel(ncores)
system.time({
  #result <- mclapply(m, hyperparameter_fx, mc.cores = ncores)
  result <- foreach(t = m, 
                    .combine = 'comb',
                    .multicombine = TRUE,
                    .init = list(list(), list(), list())) %:% 
    foreach (t2 = m, .combine = 'comb',
             .multicombine = TRUE,
             .init = list(list(), list(), list())) %dopar%{
               out <- hyperparameter_fx(t, t2)                 
             }
  
  
  
})

unlist_res <- function(res, numcol){
  unlist_res <- matrix(unlist(t(res)), ncol=numcol, byrow = TRUE)
}
gibbs_med = unlist_res(result[[1]],K)
gibbs_sd = unlist_res(result[[2]],K)
par = unlist_res(result[[3]],5)

result = list(gibbs_med=gibbs_med,gibbs_sd=gibbs_sd, par=par )
id <- which(abs(result$gibbs_med[,1] - pmin)<= 0.01 & abs(result$gibbs_med[,K] - pmax)<= 0.01)
idmax <- which.max(apply((result$gibbs_sd[id,])^2,1,sum))
hyper <- result$par[id[idmax],]

t <- hyper[1]
t2 <- hyper[2]
i <- hyper[3]
j <- hyper[4]         
l <- hyper[5]

alpha <- c(t-i,rep(l,K-2),j)
beta <- c(i,rep(t/2-l,K-2),t2-j)
