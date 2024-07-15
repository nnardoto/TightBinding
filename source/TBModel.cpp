#include "TBModel.h"

TightBinding::TightBinding()                   
{                                              
  R = arma::mat(3, 3, arma::fill::zeros);      
  K = arma::mat(3, 3, arma::fill::zeros);      
}                                              
