#include <TMB.hpp>                                // Links in the TMB libraries

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;
  
  DATA_VECTOR(X);                                     // Data vector transmitted from R
  DATA_VECTOR(cov);
  PARAMETER(beta0);                                  // Parameter value transmitted from R
  PARAMETER(beta1);
  PARAMETER(sigma);
  PARAMETER(log_kappa);
  PARAMETER_VECTOR(x);
    
  DATA_STRUCT(spde,spde_t);
  
  Type kappa = exp(log_kappa);
  
  SparseMatrix<Type> Q = Q_spde(spde, kappa);
  
  int n;  
  n = X.size();

  Type f;     // Declare the "objective function" (neg. log. likelihood)
  f= GMRF(Q)(x);
  //printf("test\n")
  //f=0;
  
  for(int i=0; i<n; i++){
    f = f -dnorm(X(i),1,sigma,true);      

  }
            // Use R-style call to normal density

  return f;
}
