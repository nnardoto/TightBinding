#include "TBModel.h"


arma::cx_mat TightBinding::GetH(double k1, double k2, double k3)
{
  // Calcula as Matrizes para um Dado ponto K (unidades de cristal)                        
  arma::cx_mat H = arma::cx_mat(nOrbitals, nOrbitals, arma::fill::zeros);                       
                                                                                           
  double h = 0.0;                                                                          
  int l, m, n;
  complex <double> Fase;
                                                                                           
  for(int i = 0; i < FockNumber; i++)                                                      
  {                                                                                        
    l = aIndex[i][0];                                                                      
    m = aIndex[i][1];                                                                      
    n = aIndex[i][2];                                                                      
                                                                                           
    // k cdot R                                                                          
    h = k1*l + k2*m + k3*n;                                                                
                                             
    Fase = cos(2.0 * PI* h) + 1i * sin(2.0 * PI* h);  
    H   += FockMatrices[i] * Fase;                                          
                                                                                           
  }                                                                                        

  H.clean(UThr);
                                                                                           
  return H;                                                                           
}

arma::vec TightBinding::BandCalc(double k1, double k2, double k3)                                                            
{                                                                                                                            
  // Calcula as Matrizes para um Dado ponto K (unidades de cristal)                                                          
  arma::cx_mat H = arma::cx_mat(nOrbitals, nOrbitals, arma::fill::zeros);                                                                         
  arma::cx_mat S = arma::cx_mat(nOrbitals, nOrbitals, arma::fill::zeros);                                                                         

  arma::vec EigVal;                                                                                                          
                                                                                                                             
                                                                                                                             
  double h = 0.0;                                                                                                            
  int l, m, n;           
  complex <double> Fase;
                                                                                                                             
  for(int i = 0; i < FockNumber; i++)
  {                                                                                                                        
    l = aIndex[i][0];                                                                                                      
    m = aIndex[i][1];                                                                                                      
    n = aIndex[i][2];                                                                                                      
                           
      
    // k cdot R at crystal cordinates
    h = k1*l + k2*m + k3*n;                                                                                                
                        
    if(OrthogonalBasis)
    {
      H += Degen[i]*FockMatrices[i] * exp(2.0i * PI * h);
    } 
    else
    {
      H += FockMatrices[i] * exp(2.0i * PI * h);
      S += Overlap[i]      * exp(2.0i * PI * h);   
    }
  }
                                                                                                                             
                                                                                                                             
  H.clean(UThr);       
  H *= 27.2;
  H.clean(UThr);

  if(OrthogonalBasis)
  {
    EigVal = arma::eig_sym(H);
  }
  else
  {
    //cout << format("{} {} {}\n", k1, k2, k3);
    S      = arma::inv(arma::sqrtmat(S));
    S.clean(UThr);
    EigVal = arma::eig_sym(S*H*S);                                                                                             
  }
  
  return EigVal;                                                                                                             
}
