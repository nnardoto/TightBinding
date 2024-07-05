#include "TBModel.h"


void TightBinding::Load(string FileName)
{
  InputDict = Parser(FileName);

  string IFormat = InputDict["IFormat"];

  if(IFormat ==    "openmx") 
  {
    cout << "==============================================\n";
    cout << "================ OPENMX MODE =================\n";
    cout << "==============================================\n\n";
    OMX_Load();
  }
  else if(IFormat == "Wannier90")
  {
    W90_Load();
  }
}


void TightBinding::MakeMap(int Neighbors)
{
  int MapSize = 2 * Neighbors + 1;                           
                                                                  
  Map = new int*[MapSize];                                        
  for(int i = 0; i < MapSize; i++)                                
    Map[i] = new int[MapSize];                                    
                                                                  
  for(int i = -Neighbors; i <= Neighbors; i++)          
  {                                                               
    for(int j = -Neighbors; j <= Neighbors; j++)        
    {                                                             
       Map[Neighbors + i][Neighbors + j] = 0;           
    }                                                             
  }                

  // Armazena internamente
  NeighborsCells = Neighbors;
}

void TightBinding::OMX_Load()
{
  // Carrega Numero de Vizinhos do Input                              
  int N = stoi(InputDict["NeighborsCells"]);        
  // Monta Mapa de Indices
  MakeMap(N);

  // Load FockMatrices                             
  OrthogonalBasis = false;
  FockCount();                                                                  
  LoadFock();                                                                   
  LoadOverlap();                                                                
}

void TightBinding::WF_SkipHead(ifstream& iFile)
{
  string Line;                                            
  // Head Reading                                         
  getline(iFile, Line);                                 
                                                          
  // Lê Vetores                                           
  for(int h = 0; h < 3; h++)  getline(iFile, Line);                            
                                                          
  iFile >> nOrbitals;                                   
  iFile >> FockNumber;                                  
                                                          
  int DD;                                     
  for(int h = 0; h < FockNumber; h++) iFile >> DD;   
}

void TightBinding::W90_Load()
{
  ifstream w90File;
  w90File.open(InputDict["HamiltonianFile"]);

  int l, m, n, i, j;
  double H1, H2;
  
  string Line;
  // Head Reading
  getline(w90File, Line);

  // Lê Vetores
  for(int h = 0; h < 3; h++)
  {
    getline(w90File, Line);
    R[h] = arma::vec(Line);
  }
 
  w90File >> nOrbitals;
  w90File >> FockNumber;

  int DD[FockNumber];
  for(int h = 0; h < FockNumber; h++) w90File >> DD[h];

  // Aloca Memoria para as Matrizes                                         
  FockMatrices     = new arma::cx_mat[FockNumber];
  arma::mat *realH = new arma::mat[FockNumber];
  arma::mat *imagH = new arma::mat[FockNumber];

        
  for(int h = 0; h < FockNumber; h++)                                       
  {                                                                         
    FockMatrices[h] = arma::cx_mat(nOrbitals, nOrbitals, arma::fill::zeros);
    realH[h]        = arma::mat(nOrbitals, nOrbitals, arma::fill::zeros);
    imagH[h]        = arma::mat(nOrbitals, nOrbitals, arma::fill::zeros);
  }                                                                         

  // Load Fock Matrices
  NeighborsCells = 8;
  MakeMap(NeighborsCells); // PRECISO DEIXAR ISSO DINAMICO
  
  int nEle = nOrbitals * nOrbitals;
  for(int h = 0; h < FockNumber; h++)
  {
    w90File >> l >> m >> n;
    Map[NeighborsCells + l][NeighborsCells + m] = 1;
    
    for(int hij = 0; hij < nEle; hij++)
    {
      w90File >> i >> j >> H1 >> H2;
    }
  }

  // Indexa indices
  int MapSize = 2 * NeighborsCells + 1;
  int mIdx = 0;
  for(int h = 0; h < MapSize; h++)
  {
    for(int g = 0; g < MapSize; g++)
    {
      if(Map[h][g])
      {
        mIdx += 1;
        Map[h][g] = mIdx;
      }
    }
  }
                                          
  aIndex = new int*[FockNumber];        
  for(int h = 0; h < FockNumber; h++)   
  {                                     
    aIndex[h] = new int[3];             
  }                                     

  // Finalmente Lê a Hamiltoniana
  w90File.close();                               
  w90File.open(InputDict["HamiltonianFile"]);    
  WF_SkipHead(w90File);                                 

  for(int h = 0; h < FockNumber; h++)                      
  {                                                        
    w90File >> l >> m >> n;                                
                                                           
    for(int hij = 0; hij < nEle; hij++)                    
    {                                                      
      w90File >> i >> j >> H1 >> H2;
      
      aIndex[h][0] = l;
      aIndex[h][1] = m;
      aIndex[h][2] = n;
      
      realH[Index(l, m, n)](i - 1, j - 1) = H1;
      imagH[Index(l, m, n)](i - 1, j - 1) = H2;
    }                                                      
  }                                                        

  for(int h = 0; h < FockNumber; h++)
  {
    FockMatrices[h] = arma::cx_mat(realH[h], imagH[h]);
  }
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
  FockMatrices = new arma::cx_mat[FockNumber];                                                 
                        
  cout << "========= ALOCATION OF FOCK MATRICES =========\n";
  for(int i = 0; i < FockNumber; i++)                                                       
  {                             
    cout << format("Alocation of {: 4d}/{: 4d} Fock Matrices", i + 1, FockNumber);
    FockMatrices[i] = arma::cx_mat(nOrbitals, nOrbitals, arma::fill::zeros);                          
    cout << format(" ... ok\n");
  }                                                                                
  cout << "==============================================\n\n";
                                                                                   
  ifstream Hamiltonian;
  Hamiltonian.open(InputDict["HamiltonianFile"]);

  int l, m, n, i, j;
  double U;
  
  cout << "============ LOAD OF FOCK MATRICES ===========\n";
  while(Hamiltonian >> l >> m >> n >> i >> j >> U)
  {
    FockMatrices[Index(l, m, n)](i - 1, j - 1) = U;
  }
  cout << format("{: 4d} Matrices Loaded\n", FockNumber);
  cout << "==============================================\n\n";

  Hamiltonian.close();
}                                                                                  
                                                                                   

void TightBinding::LoadOverlap()                   
{
  // Aloca Memoria para as Matrizes                                                      
  Overlap = new arma::cx_mat[FockNumber];                                              
                                                                                         
  cout << "======== ALOCATION OF OVERLAP MATRICES =======\n";
  for(int i = 0; i < FockNumber; i++)                                                    
  {                                                                                      
    cout << format("Alocation of {: 4d}/{: 4d} Overlap Matrices", i + 1, FockNumber);
    Overlap[i] = arma::cx_mat(nOrbitals, nOrbitals, arma::fill::zeros);
    cout << format(" ... ok\n");
  }                                                         
  cout << "==============================================\n\n";
                                                                                         
  ifstream OLP;                                                                  
  OLP.open(InputDict["OverlapFile"]);                                        
                                                                                         
  int l, m, n, i, j;                                                                     
  double U;                                                                              
                                      

  cout << "========== LOAD OF OVERLAP MATRICES ==========\n";
  while(OLP >> l >> m >> n >> i >> j >> U)                                       
  {                                                                                      
    Overlap[Index(l, m, n)](i - 1, j - 1) = U;                                      
  }                                                                                      
  cout << format("{: 4d} Matrices Loaded\n", FockNumber);
  cout << "==============================================\n\n";
  
  OLP.close();                                                                   
}
