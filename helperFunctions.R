DLT_status <- function(data1,assigned_dose_idx,cohort){
  if (cohort==1) {
    num = 2
  }else if (cohort==2){
    num = 2
  }else{
    num =1
    }
  z1 <- sum(rbinom(num,1,data1$p[assigned_dose_idx[1]]))
  z2 <- sum(rbinom(num,1,data1$p[assigned_dose_idx[2]]))
  dlt<-c(z1,z2)
  return(dlt)
}


AssignDose_basedDirection <-function(assigned_dose,cohort,x1,y1,posterior,position_in_cohort){
  xh = assigned_dose[[cohort-1]][[position_in_cohort]][1]
  Dh = posterior[xh,] 
  yh = y1[which.min(abs(Dh-theta))]
  if (yh-assigned_dose[[cohort-1]][[position_in_cohort]][2] >1) yh = assigned_dose[[cohort-1]][[position_in_cohort]][2]+1  #no skipping dose while escalation
  
  yv = assigned_dose[[cohort-1]][[position_in_cohort]][2]
  Dv = posterior[,yv] 
  xv = x1[which.min(abs(Dv-theta))]
  if (xv- assigned_dose[[cohort-1]][[position_in_cohort]][1] >1) xv = assigned_dose[[cohort-1]][[position_in_cohort]][1]+1  #no skipping dose while escalation
  direction = c("horiz","vert")
  
  if (min(Dh) > 1.5*theta){
    xh = xv
    yh = yv
   }else if(min(Dv) > 1.5*theta){
    xv = xh
    yv = yh
   }
  
  if(min(Dh) > 1.5*theta & min(Dv) > 1.5*theta){
    id.min <- which.min(c(min(Dh),min(Dv)))
    if (id.min==1){
      xh = xv
      yh = yv
    }else{
      xv = xh
      yv = yh
    }
  }
  
  if (cohort==2 || all(unlist(assigned_dose[[cohort-1]])==1)){
    return(list(c(xh,yh),c(xv,yv)))
  }else{

      frac_h = frac_v = 1/2
    }
    
    direct = sample(direction,1,replace=TRUE, prob = c(frac_v,frac_h))
    
    if (direct=="horiz"){
      return(c(xh,yh))
    }else{
      return(c(xv,yv))
    }

}
