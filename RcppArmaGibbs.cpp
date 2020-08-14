#include <RcppArmadillo.h>
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
// double rTruncBetaC3(double alpha, double beta, double lbound, double ubound) {
//   double F0 = R::pbeta(lbound, alpha, beta,1,1);
//   double F1 = R::pbeta(ubound, alpha, beta,1,1);
//   double U = R::runif(max(F0,-5000.00),F1);
//   double r = R::qbeta(U, alpha, beta,1,1);
//   return r;
// }

double rTruncBetaC2(double alpha, double beta, double lbound, double ubound) {
  double S0 = R::pbeta(lbound, alpha, beta, 0,0);
  double S1 = R::pbeta(ubound, alpha, beta, 0,0);
  double U = R::runif(S1, S0);
  double r = R::qbeta(U, beta, alpha,0,0);
  return r;
}
// [[Rcpp::export]]
double rTruncBetaC(double alpha, double beta, double lbound, double ubound) {
  double F0 = R::pbeta(lbound, alpha, beta,1,0);
  double F1 = R::pbeta(ubound, alpha, beta,1,0);
  double U = R::runif(F0,F1);
  //double r = (max(U)==1)*rTruncBetaC2(alpha, beta, lbound, ubound) + R::qbeta(U, alpha, beta,1,0);
  double r = R::qbeta(U, alpha, beta,1,0);
  return r;
}


// [[Rcpp::export]]
mat latticeGibbsSamplerC(mat A, vec alpha, vec beta, int nSamp, int burnT){
  int K = alpha.n_elem;
  
  //int K = alpha.size()
  //int N1= A.n_rows, N2 = A.n_cols;
  //int K = N1*N2;
  umat edgelist = ind2sub(size(A),find(A!=0));
  urowvec r = edgelist.row(0);
  urowvec c = edgelist.row(1);
  int Tmax = nSamp + burnT;
  mat lattice_sample = zeros(nSamp, K);
  mat bounds = zeros(K,2);
  rowvec cur_sample(K);
  cur_sample.fill(1./K);
  
  for (int t = 1; t<=Tmax; t++) {
    for (int i=0; i< K; i++){
      uvec idx_r = find(r==i);  
      if(!idx_r.is_empty()){
        uvec vecc_out = c(idx_r);
        bounds(i,1) = min(cur_sample(vecc_out));
      }else{
        bounds(i,1) = 1;
      }
      uvec idx_c = find(c==i);  
      if(!idx_c.is_empty()){
        uvec vecc_out = r(idx_c);
        bounds(i,0) = max(cur_sample(vecc_out));
      }else{
        bounds(i,0) = 0;
      }
      cur_sample(i) = rTruncBetaC(alpha(i), beta(i), bounds(i,0), bounds(i,1));
      while((cur_sample(i) < bounds(i,0)) | (cur_sample(i) >bounds(i,1))){
        cur_sample(i) = rTruncBetaC(alpha(i), beta(i), bounds(i,0), bounds(i,1));
      }
    }
    if (t > burnT) {
      lattice_sample.row(t-burnT-1) = cur_sample;
    }
  }
  return(lattice_sample);
}



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//


