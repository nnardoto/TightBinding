#include "BandCalc.h"

void TightBinding::SetR(string FileName)
{
  arma::mat RR;
  RR.load(FileName);

  R[0] = RR.row(0).t();
  R[1] = RR.row(1).t();
  R[2] = RR.row(2).t();


  // Calculate Reciprocal Space Vectors
  double V = abs(arma::dot(arma::cross(R[0], R[1]), R[2])); 
  
  K[0] = 2*PI*arma::cross(R[1], R[2]) / V;
  K[1] = 2*PI*arma::cross(R[2], R[0]) / V;
  K[2] = 2*PI*arma::cross(R[0], R[1]) / V;
}


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

int TightBinding::LoadFockMatrices(int N)
{
  char iName[80];
  int Index = 0;
  FockNumber = N;

  // Aloca Memoria para as Matrizes
  FockMatrices = new arma::mat[9];

  for(int i = 0; i < 9; i++)
  {
    FockMatrices[i] = arma::mat(8, 8, arma::fill::zeros);
  }


  for(int i = -1; i <= 1; i++)
  {
    for (int j = -1; j <= 1; j++)
    {
      sprintf(iName, "H%i%i0.dat", i, j);
      FockMatrices[Index].load(iName);
      Rn[Index][0] = i;
      Rn[Index][1] = j;
      Rn[Index][2] = 0;
      Index += 1;
      // LEMBRAR DE LIMPAR A MEMORIA NA SAIDA
    }
  }

  return Index;
}

                                                            
int TightBinding::LoadOverlap(int N)                   
{                                                           
  char iName[80];                                           
  int Index = 0;                                            
  FockNumber = N;                                           
                                                            
  // Aloca Memoria para as Matrizes                         
  Overlap = new arma::mat[9];                         
                                                            
  for(int i = 0; i < 9; i++)                               
  {                                                         
    Overlap[i] = arma::mat(8, 8, arma::fill::zeros);   
  }                                                         
                                                            
                                                            
  for(int i = -1; i <= 1; i++)                              
  {                                                         
    for (int j = -1; j <= 1; j++)                           
    {                                                       
      sprintf(iName, "S%i%i0.dat", i, j);                   
      Overlap[Index].load(iName);                      
      Rn[Index][0] = i;                                     
      Rn[Index][1] = j;                                     
      Rn[Index][2] = 0;                                     
      Index += 1;                                           
      // LEMBRAR DE LIMPAR A MEMORIA NA SAIDA               
    }                                                       
  }                                                         
                                                            
  return Index;                                             
}
