#include <TMB.hpp>                                // Links in the TMB libraries

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(x);                                     // Data vector transmitted from R
  DATA_VECTOR(cov);
  PARAMETER(beta0);                                  // Parameter value transmitted from R
  PARAMETER(beta1);
  PARAMETER(sigma)
  
  int n;  
  n = x.size();

  Type f;     // Declare the "objective function" (neg. log. likelihood)
  f=0;
  for(int i=0; i<n; i++){
    f = f -dnorm(x[i],beta0 + beta1*cov[i],sigma,true);       
  }
            // Use R-style call to normal density

  return f;
}
