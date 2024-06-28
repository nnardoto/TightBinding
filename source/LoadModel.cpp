#include "TBModel.h"


void TightBinding::Load(string FileName)
{
  InputDict = Parser(FileName);

  // Valida minimanente o input
  NeighborsCells = stoi(InputDict["NeighborsCells"]);
 
  int MapSize = 2 * NeighborsCells + 1;
  
  Map = new int*[MapSize];
  for(int i = 0; i < MapSize; i++)
    Map[i] = new int[MapSize];

  for(int i = -NeighborsCells; i <= NeighborsCells; i++)
  {
    for(int j = -NeighborsCells; j <= NeighborsCells; j++)
    {
        Map[NeighborsCells + i][NeighborsCells + j] = 0;
    }
  }

  // Load FockMatrices  
  FockCount();
  LoadFock();
  LoadOverlap();
}

int TightBinding::Index(int i, int j, int k)
{
  int ii = NeighborsCells + i;
  int jj = NeighborsCells + j;
  int kk = NeighborsCells + k;
  
  return Map[ii][jj] - 1;
}


void TightBinding::FockCount()
{
  ifstream Hamiltonian;                                  
  Hamiltonian.open(InputDict["HamiltonianFile"]);        
                                                         
  int l, m, n, i, j;                                     
  double U;                                              
                                                         
  while(Hamiltonian >> l >> m >> n >> i >> j >> U)       
  {                                                      
    Map[NeighborsCells + l][NeighborsCells + m] = 1;      
  }                                          

  int MapSize = 2 * NeighborsCells + 1;
  for(int i = 0; i < MapSize; i++)
    for(int j = 0; j < MapSize; j++)
    {
      if(Map[i][j])
      {
        FockNumber += 1;
        Map[i][j] = FockNumber;
      }
    }
  
    aIndex = new int*[FockNumber];                          
    for(int i = 0; i < FockNumber; i++)                     
    {                                                       
      aIndex[i] = new int[3];                               
    }                         

  Hamiltonian.clear();
  Hamiltonian.seekg(0);

  while(Hamiltonian >> l >> m >> n >> i >> j >> U)     
  {                                                    
    aIndex[Index(l,m,n)][0] = l;    
    aIndex[Index(l,m,n)][1] = m;    
    aIndex[Index(l,m,n)][2] = n;    
  }                                                    

  Hamiltonian.close();  
}

void TightBinding::LoadFock()                                         
{ 
  // Aloca Memoria para as Matrizes
  FockMatrices = new arma::mat[FockNumber];                                                 
                               
  for(int i = 0; i < FockNumber; i++)                                                       
  {                             
    FockMatrices[i] = arma::mat(8, 8, arma::fill::zeros);                          
  }                                                                                
                                                                                   
  ifstream Hamiltonian;
  Hamiltonian.open(InputDict["HamiltonianFile"]);

  int l, m, n, i, j;
  double U;
  
  while(Hamiltonian >> l >> m >> n >> i >> j >> U)
  {
    FockMatrices[Index(l, m, n)](i - 1, j - 1) = U;
  }

  Hamiltonian.close();
}                                                                                  
                                                                                   

void TightBinding::LoadOverlap()                   
{
  // Aloca Memoria para as Matrizes                                                      
  Overlap = new arma::mat[FockNumber];                                              
                                                                                         
  for(int i = 0; i < FockNumber; i++)                                                    
  {                                                                                      
    Overlap[i] = arma::mat(8, 8, arma::fill::zeros);                                
  }                                                                                      
                                                                                         
  ifstream OLP;                                                                  
  OLP.open(InputDict["OverlapFile"]);                                        
                                                                                         
  int l, m, n, i, j;                                                                     
  double U;                                                                              
                                                                                         
  while(OLP >> l >> m >> n >> i >> j >> U)                                       
  {                                                                                      
    Overlap[Index(l, m, n)](i - 1, j - 1) = U;                                      
  }                                                                                      
                                                                                         
  OLP.close();                                                                   
}
