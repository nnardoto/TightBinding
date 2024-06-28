#include "TBModel.h"

arma::vec TightBinding::BandCalc(double k1, double k2, double k3)                                                            
{                                                                                                                            
  // Calcula as Matrizes para um Dado ponto K (unidades de cristal)                                                          
  arma::mat H1 = arma::mat(8, 8, arma::fill::zeros);                                                                         
  arma::mat H2 = arma::mat(8, 8, arma::fill::zeros);                                                                         
  arma::mat S1 = arma::mat(8, 8, arma::fill::zeros);                                                                         
  arma::mat S2 = arma::mat(8, 8, arma::fill::zeros);                                                                         
                                                                                                                             
  arma::vec EigVal;                                                                                                          
                                                                                                                             
  arma::cx_mat H, S;                                                                                                         
                                                                                                                             
  arma::vec k = arma::vec();                                                                                                 
  arma::vec r = arma::vec();                                                                                                 
                                                                                                                             
                                                                                                                             
  double h = 0.0;                                                                                                            
  int l, m, n;                                                                                                               
                                                                                                                             
  int Index = 0;                                                                                                             
  for(int i = -1; i <= 1; i++)                                                                                               
  {                                                                                                                          
    for(int j = -1; j <=1; j++)                                                                                              
    {                                                                                                                        
      l = Rn[Index][0];                                                                                                      
      m = Rn[Index][1];                                                                                                      
      n = Rn[Index][2];                                                                                                      
                                                                                                                             
      h = k1*l + k2*m + k3*n;                                                                                                
                                                                                                                             
      H1 += FockMatrices[Index] * cos(2.0 * PI* h);                                                                          
      H2 += FockMatrices[Index] * sin(2.0 * PI* h);                                                                          
      Index += 1;                                                                                                            
    }                                                                                                                        
  }                                                                                                                          
                                                                                                                             
  Index = 0;                                                                                                                 
  for(int i = -1; i <= 1; i++)                                                                                               
  {                                                                                                                          
    for(int j = -1; j <= 1; j++)                                                                                             
    {                                                                                                                        
      l = Rn[Index][0];                                                                                                      
      m = Rn[Index][1];                                                                                                      
      n = Rn[Index][2];                                                                                                      
                                                                                                                             
      h = k1*l + k2*m + k3*n;                                                                                                
                                                                                                                             
      S1 += Overlap[Index] * cos(2.0*PI*h);                                                                                  
      S2 += Overlap[Index] * sin(2.0*PI*h);                                                                                  
      Index += 1;                                                                                                            
    }                                                                                                                        
  }                                                                                                                          
                                                                                                                             
  H = arma::cx_mat(H1, H2);                                                                                                  
  S = arma::cx_mat(S1, S2);                                                                                                  
                                                                                                                             
  S = inv(sqrtmat(S));                                                                                                       
                                                                                                                             
  EigVal = arma::eig_sym(S*H*S);                                                                                             
                                                                                                                             
  return EigVal;                                                                                                             
}                                                                                                                            
