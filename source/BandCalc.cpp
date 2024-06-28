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
                                                                                                                             
  double h = 0.0;                                                                                                            
  int l, m, n;                                                                                                               
                                                                                                                             
  for(int i = 0; i < FockNumber; i++)
    {                                                                                                                        
      l = aIndex[i][0];                                                                                                      
      m = aIndex[i][1];                                                                                                      
      n = aIndex[i][2];                                                                                                      
                                                  
      // k cdot R
      h = k1*l + k2*m + k3*n;                                                                                                
                                                                                                                             
      H1 += FockMatrices[i] * cos(2.0 * PI* h);                                                                          
      H2 += FockMatrices[i] * sin(2.0 * PI* h);                                                                          
    }                                                                                                                        
                                                                                                                             
    for(int i = -1; i <= FockNumber; i++)                                                                                             
    {                                                                                                                        
      l = aIndex[i][0];                                                                                                      
      m = aIndex[i][1];                                                                                                      
      n = aIndex[i][2];                                                                                                      
                                                                                                                             
      h = k1*l + k2*m + k3*n;                                                                                                
                                                                                                                             
      S1 += Overlap[i] * cos(2.0*PI*h);                                                                                  
      S2 += Overlap[i] * sin(2.0*PI*h);                                                                                  
    }                                                                                                                        
                                                                                                                             
  H = arma::cx_mat(H1, H2);                                                                                                  
  S = arma::cx_mat(S1, S2);                                                                                                  
                                                                                                                             
  S = inv(sqrtmat(S));                                                                                                       
                                                                                                                             
  EigVal = arma::eig_sym(S*H*S);                                                                                             
                                                                                                                             
  return EigVal;                                                                                                             
}
                                                                                 
arma::mat TightBinding::PathCalc(arma::vec Origin, arma::vec Destination, int N) 
{                                                                                
  // Armazena estrutura de bandas em uma matriz                                  
  arma::mat BandStructure(nOrbitals, N, arma::fill::zeros);                      
  arma::vec step = (Destination - Origin)/N;                                     
  arma::vec kp = Origin;                                                         
                                                                                 
  for(int i = 0; i <= N; i++)                                                    
  {                                                                              
    BandStructure.col(i) = BandCalc(0.0, 0.0, 0.0);                                  
    kp += step;                                                                  
  }                                                                              
                                                                                 
  return BandStructure;                                                          
}                                                                                

