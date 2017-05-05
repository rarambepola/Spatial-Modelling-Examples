#include <TMB.hpp>                                // Links in the TMB libraries

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;
  
  DATA_VECTOR(X);                                     // Data vector transmitted from R
  DATA_VECTOR(cov);

  PARAMETER(beta5);                                  // Parameter value transmitted from R
  PARAMETER(beta1);
  PARAMETER(sigma);
  PARAMETER(log_kappa);
  PARAMETER_VECTOR(x);
  
  DATA_SPARSE_MATRIX(Apixel);
  DATA_IVECTOR(box_total);
    
  DATA_STRUCT(spde,spde_t);
  
  Type kappa = exp(log_kappa);
  
  Type log_kappa_mean = -3.0;
  Type log_kappa_sd = 1.0;
  Type beta_mean = 0.0;
  Type beta_sd = 1.0;
  
  SparseMatrix<Type> Q = Q_spde(spde, kappa);
  
  int n;  
  n = X.size();
  
  Type f;     // Declare the "objective function" (neg. log. likelihood)
  f= GMRF(Q)(x);
  //printf("test\n")
  //f=0;
  
  f = f - dnorm(beta5, beta_mean, beta_sd);
  f = f - dnorm(beta1, beta_mean, beta_sd);
  f = f - dnorm(log_kappa, log_kappa_mean, log_kappa_sd);
  
  vector<Type> field = Apixel*x;
  
  int mstart = 0;
  int mstop = 0;
  
  for(int i=0; i<n; i++){
    Type total=0;
    mstop = mstart + box_total(i);
    for(int j=mstart; j<mstop; j++){
      total = total + beta5 + beta1*cov(j) + field(j);
    }
    f = f -dnorm(X(i), total,sigma,true); 
    mstart = mstop;
  }
            // Use R-style call to normal density

  return f;
}
