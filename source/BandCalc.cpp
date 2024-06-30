#include "TBModel.h"


arma::cx_mat TightBinding::GetH(double k1, double k2, double k3)
{
  // Calcula as Matrizes para um Dado ponto K (unidades de cristal)                        
  arma::mat H1 = arma::mat(nOrbitals, nOrbitals, arma::fill::zeros);                       
  arma::mat H2 = arma::mat(nOrbitals, nOrbitals, arma::fill::zeros);                       
                                                                                           
  arma::cx_mat H;                                                                    
                                                                                           
  double h = 0.0;                                                                          
  int l, m, n;                                                                             
                                                                                           
  for(int i = 0; i < FockNumber; i++)                                                      
  {                                                                                        
    l = aIndex[i][0];                                                                      
    m = aIndex[i][1];                                                                      
    n = aIndex[i][2];                                                                      
                                                                                           
//      cout << l << ' ' << m << ' ' << n << endl;                                         
                                                                                           
      // k cdot R                                                                          
    h = k1*l + k2*m + k3*n;                                                                
                                                                                           
    H1 += FockMatrices[i] * cos(2.0 * PI* h);                                              
    H2 += FockMatrices[i] * sin(2.0 * PI* h);                                              
                                                                                           
  }                                                                                        
                           

  H = arma::cx_mat(H1, H2);
  H.clean(UThr);
                                                                                           
                                                                                           
  return H;                                                                           
}

arma::vec TightBinding::BandCalc(double k1, double k2, double k3)                                                            
{                                                                                                                            
  // Calcula as Matrizes para um Dado ponto K (unidades de cristal)                                                          
  arma::mat H1 = arma::mat(nOrbitals, nOrbitals, arma::fill::zeros);                                                                         
  arma::mat H2 = arma::mat(nOrbitals, nOrbitals, arma::fill::zeros);                                                                         
  arma::mat S1 = arma::mat(nOrbitals, nOrbitals, arma::fill::zeros);                                                                         
  arma::mat S2 = arma::mat(nOrbitals, nOrbitals, arma::fill::zeros);                                                                         
                                                                                                                             
  arma::vec EigVal;                                                                                                          
                                                                                                                             
  arma::cx_mat H, S;                                                                                                         
                                                                                                                             
  double h = 0.0;                                                                                                            
  int l, m, n;                                                                                                               
                                                                                                                             
  for(int i = 0; i < FockNumber; i++)
  {                                                                                                                        
    l = aIndex[i][0];                                                                                                      
    m = aIndex[i][1];                                                                                                      
    n = aIndex[i][2];                                                                                                      
                           
//      cout << l << ' ' << m << ' ' << n << endl;
      
      // k cdot R
    h = k1*l + k2*m + k3*n;                                                                                                
                                                                                                                             
    H1 += FockMatrices[i] * cos(2.0 * PI* h);                                                                          
    H2 += FockMatrices[i] * sin(2.0 * PI* h);                                                                          
      
    S1 += Overlap[i] * cos(2.0*PI*h);   
    S2 += Overlap[i] * sin(2.0*PI*h);   
  }                                                                                                                        
                                                                                                                             
                                                                                                                             
  H = arma::cx_mat(H1, H2).clean(UThr);                                                                                                  
  S = arma::cx_mat(S1, S2).clean(UThr);                                                                                                  

  S = inv(sqrtmat(S));                                                                                                       
                                                                                                                             
  EigVal = arma::eig_sym(S*H*S);                                                                                             
                                                                                                                             
  return EigVal;                                                                                                             
}
                                                                                 
void TightBinding::PathCalc() 
{                                                                                
  // Lê Input para Montar Caminho
  ifstream fp;
  int nPath = stoi(InputDict["nPath"]);
  
  // Armazena Caminho
  arma::vec k[nPath], kk[nPath];
  for(int i = 0; i < nPath; i++)
  {
    k[i]  = arma::vec(3, arma::fill::zeros);
    kk[i] = arma::vec(3, arma::fill::zeros);
  }

  // abre arquivo
  fp.open(InputDict["KPathFile"]);

  string HighA[nPath], HighB[nPath];
  int N[nPath];
  for(int i = 0; i < nPath; i++)
  {
    fp >> HighA[i] >>  k[i](0) >>  k[i](1) >>  k[i](2) 
       >> HighB[i] >> kk[i](0) >> kk[i](1) >> kk[i](2) >> N[i];
  }

  fp.close();
  // Armazena estrutura de bandas em uma matriz
  int kpoints = 0;
  for(int i = 0; i < nPath; i++)
  {
    kpoints += N[i];
  }

  arma::mat BandStructure(nOrbitals, kpoints, arma::fill::zeros);
  arma::vec FullPath(kpoints, arma::fill::zeros);
  arma::vec kp_ant(3, arma::fill::zeros);

  int column = 0;
  double kkp = 0.0;
  kp_ant = k[0];
  for(int i = 0; i < nPath; i++)
  {
    arma::vec step = (kk[i] - k[i])/(N[i]);                    
    arma::vec kp = kp_ant;                                           

   // cout << "========================\n";                     
   // cout << "Start: " << k[i](0) << ' ' << k[i](1) << ' ' << k[i](2) << endl;
   // cout << "kp: " << kp(0) << ' ' << kp(1) << ' ' << kp(2) << endl;
   // cout << "step: " << step(0) << ' ' << step(1) << ' ' << step(2) << endl;
   // cout << "End: " << kk[i](0) << ' ' << kk[i](1) << ' ' << kk[i](2) << endl;
   // cout << "========================\n";                     

    for(int j = 0; j < N[i]; j++)                                                    
    {        
      BandStructure.col(column) = BandCalc(kp(0), kp(1), kp(2));                                  
      kkp += norm(step);
      FullPath(column) = kkp;
      kp += step; 
     // cout << kp.t();
      column += 1;
    }
    kp_ant = kk[i];
  }

  // Imprime Estrutura de Bandas
  for(int i = 0; i < nOrbitals; i++)
  {
    for(int j = 0; j < kpoints; j++)
    {
      cout << FullPath(j) << "\t\t" << BandStructure.col(j)(i) << endl;
    }
    cout << endl;
  }
}                                                                                

