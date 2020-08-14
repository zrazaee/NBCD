source("helperFunctions.R")
source('GibbsPosterior.R')
require(parallel)
NBCD  <- function(d,N,alpha,beta, theta, Dlt.true, ntrial = 2000, gamma = 0.1, eps = 0.8, nSamp =  10000, burnT = 1000, ncores = 4){
    
        K <- prod(d)
       
        ncohort <- floor((N-8)/2+2)

        A <- matrix(0,K,K)
        for (i in 1:(K-1)){
            if  (i%%d[2]!=0) {A[i,i+1]= 1}
            if  (i<= (d[2]*(d[1]-1))){
                A[i,i+d[2]]=1
            }
        }
        
       
       Rcpp::sourceCpp('RcppArmaGibbs.cpp')
       

       num = sum(alpha)+sum(beta)
       
       x1 = seq(1,d[1],length.out=d[1])
       y1 = seq(1,d[2],length.out=d[2])
       
       dose_grid <- data.frame(x=rep(x1,each=d[2]),y=rep(y1,d[1]))
       data1 = data.frame(dose_grid,true_posterior=NA)
       data1$p =  as.vector(t(Dlt.true))
       
       
       DLT <- rep(0,N)   # DLT status for each patient
       
       result = list()
       num_vec = rep(0,ntrial)
       terminated = 0
       terminated_id = c()
       assigned_dose = list()
       assigned_dose_mat = DLT_mat = mat.or.vec(ncohort,2)
       z = nn = rep(0,K)
       patient_dose <- mat.or.vec(N,K)  # which patient is assigned to which dose
       pat_dose_DLT <- data.frame(cbind(patient_dose,DLT))
       
       trial_fx <- function(trial){
           for (cohort in 1:ncohort){
               
               if (cohort==1){
                   assigned_dose_idx = c(1,1)
                   p1_dose = p2_dose = c(1,1)
                   assigned_dose[[cohort]] = list(p1_dose, p2_dose)
                   
               }else if (cohort==2 || sum(assigned_dose_idx)==2){
                   
                   temp = AssignDose_basedDirection(assigned_dose,cohort,x1,y1,posterior,1)
                   p1_dose = temp[[1]]
                   p2_dose = temp[[2]]
                   
               }else{
                   
                   p1_dose = AssignDose_basedDirection(assigned_dose,cohort,x1,y1,posterior,1)
                   p2_dose = AssignDose_basedDirection(assigned_dose,cohort,x1,y1,posterior,2)
                   
                   
               }
               p1_idx = (p1_dose[1]-1)*d[2] + p1_dose[2]
               p2_idx = (p2_dose[1]-1)*d[2] + p2_dose[2]
               
               assigned_dose[[cohort]] = list(p1_dose, p2_dose )
               assigned_dose_idx = c(p1_idx,p2_idx)
               assigned_dose_mat[cohort,] = assigned_dose_idx
               
               pat_dose_DLT[2*cohort-1,assigned_dose_idx[1]] = (cohort==1)*1+(cohort==2)*1+1
               pat_dose_DLT[2*cohort,assigned_dose_idx[2]] = (cohort==1)*1+(cohort==2)*1+1
          
               dlt = DLT_status(data1,assigned_dose_idx,cohort)
               pat_dose_DLT[2*cohort-1,"DLT"] = dlt[1]
               pat_dose_DLT[2*cohort,"DLT"] = dlt[2]
               nn = colSums(pat_dose_DLT[,(1:K)])
               z[assigned_dose_idx[1]] = z[assigned_dose_idx[1]]+ dlt[1]
               z[assigned_dose_idx[2]] = z[assigned_dose_idx[2]]+ dlt[2]
               DLT_mat[cohort,] = dlt

               w = 1+2*num/sum(nn)
               
               gibbs_sample = GibbsPosterior(d,alpha,beta,w*nn,w*z,nSamp, burnT)
               posterior = matrix(gibbs_sample$median,d[1],d[2],byrow=TRUE)
               
               if ((sum(nn) >= .1*N) & mean(gibbs_sample$lattice_sample[,1] > (theta + gamma), na.rm=TRUE) > eps){
                   terminated_id = rbind(terminated_id, trial)
                   terminated = terminated +1
                   break}
               
           }
           res = list(DLT_mat=DLT_mat, assigned_dose_mat = assigned_dose_mat, posterior=posterior, z=z, nn=nn, terminated_id = terminated_id)
           
       }
       system.time({
           result <- mclapply(1:ntrial, trial_fx, mc.cores = ncores)
       })
     return(result)
    }
################################################

d <- c(4,4)

t <- 5.2653846
t2 <- 13.9730769
i <- 0.7405983
j <- 0.2
l <- 0.4

K <-  prod(d)
alpha <- c(t-i,rep(l,K-2),j)
beta <- c(i,rep(t/2-l,K-2),t2-j)

theta <- 0.2
N <- 50
ntrial <- 2
delta <- 0.1

Dlt.true <- 0.01 * matrix(c( 1, 2,  3,  4,
                            4, 10, 15, 20,
                            6, 15, 30, 45,
                            10, 30, 50, 80),nrow=4,ncol=4,byrow=T)

result = NBCD(d,N,alpha,beta, theta, Dlt.true, ntrial)


trialmtd <-  mtd_index <- list()
patientcounts <- patienttox <- 0

termin = c()
for (i in 1:ntrial){
  termin = c(termin, result[[i]]$terminated_id)
}

for (m in 1:ntrial) {
    patientcounts = patientcounts + result[[m]]$nn
    patienttox = patienttox + result[[m]]$z
    trialmtd[[m]] <- result[[m]]$posterior
    flag = 0
    
    if (any(termin%in%m)){
        mtd_index[[m]] = 0
    }else{
        x = as.vector(t(trialmtd[[m]]))-theta
        if(mean(trialmtd[[m]]>theta)>=K/2) {flag=1}
        l = 0.05
        u = 0
        temp = 0
        
        while(all(temp==0) &((l<= 0.1) ||(u<= 0.05) )){
            temp = (x >=-l & x<=u)*(result[[m]]$nn)*1
            l = l+ (l<= 0.1)*0.05
            u = u + (u<= 0.05)*(0.01*flag +0.025*(1-flag))
            
        }
        
        if(all(temp==0)){
            mtd_index[m]=0
        }else{
            
            mtd_index[[m]]=which(temp>1)
            
            if(all(temp<=1)){
                temp = (x >=-0.1 & x<=0)*(result[[m]]$nn)*1
                mtd_index[[m]]=which(temp==1)
            }
        }
        
    }
}
tab= table(unlist(mtd_index))
t = tab/sum(tab)*100

Gtrue<- ( abs(Dlt.true - theta) <= delta)
gt = as.vector(t(Gtrue))
s = unique(unlist((mtd_index)))
s = s[s>0]
gts = gt[s]*s
t[as.character(gts[gts>0])]
MTD_select <- sum(t[as.character(gts[gts>0])])  # MTD selection percentage

gttrue = as.vector(t(Gtrue))
alloc_percent <- sum( 100*(patientcounts*gttrue)/(N*ntrial) )   # allocation percentage within 10% of theta







