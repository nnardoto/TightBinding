#include "TBModel.h"

                                                                            
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

  cout << "==============================================\n\n";
  cout << "============= KPATH CALCULATION ==============\n";  
  cout << "==============================================\n\n";
  int column = 0;                                                           
  double kkp = 0.0;                                                         
  kp_ant = k[0];                                                            
  for(int i = 0; i < nPath; i++)                                            
  {                                                                         
    arma::vec step = (kk[i] - k[i])/(N[i]);                                 
    arma::vec kp = kp_ant;                                                  
                                                                            
    for(int j = 0; j < N[i]; j++)                                           
    { 
      BandStructure.col(column) = BandCalc(kp(0), kp(1), kp(2));         
      kkp += norm(step);                                                    
      FullPath(column) = kkp;                                               
      kp += step;                                                           
      column += 1;                                                          
      cout << j << ' ';
    }                                                                       
    cout << endl;
    kp_ant = kk[i];                                                         
  }
  cout << "==============================================\n\n";

  // Imprime Estrutura de Bandas                                           
  ofstream OutBand;
  OutBand.open("Bands.dat");

  for(int i = 0; i < nOrbitals; i++)                                        
  {                                                                         
    for(int j = 0; j < kpoints; j++)                                        
    {                                                                       
      OutBand << format("{: 6.4f}\t{: 8.4f}\n", FullPath(j), BandStructure.col(j)(i));
    }                                                                       
      OutBand << endl;                                                           
  }                        
  OutBand.close();
}                                                                           
                                                                            
