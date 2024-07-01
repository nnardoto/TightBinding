#include "TBModel.h"
#include <format>
#define TOEV  27.2114079527

void TightBinding::Export2Xatu()
{
  int Dimensions = stoi(InputDict["Dimensions"]);
  // int nOrbitals  = GetOrbitals();
  arma::mat R, Atoms;
  arma::vec newR;
  R.load(InputDict["CellFile"]);
  Atoms.load(InputDict["AtomsFile"]);
  
  cout << "# dimension" << endl;
  cout << Dimensions << endl;

  cout << "# norbitals" << endl;
  cout << "4 4" << endl;
  
  cout << "# bravaislattice" << endl;
  for(int i = 0; i < Dimensions; i++)
  {
    cout << std::format("{: 11.8f}    {: 11.8f}    {: 11.8f}\n", R.row(i)(0), R.row(i)(1), R.row(1)(2));
  }

  cout << "# motif" << endl;
  for(int i = 0; i < 2; i++)
  {
    cout << std::format("{: 11.8f}    {: 11.8f}    {: 11.8f}    {: 11.8f}\n", Atoms.row(i)(0), Atoms.row(i)(1), Atoms.row(i)(2), Atoms.row(i)(3));
  }

  cout << "# bravaisvectors" << endl;
  for(int i = 0; i < FockNumber; i++)
  {
    arma::vec A(3);
    A(0) = aIndex[i][0];
    A(1) = aIndex[i][1];
    A(2) = aIndex[i][2];
    newR = R * A;
    cout << std::format("{: 11.8f}    {: 11.8f}    {: 11.8f}\n", newR(0), newR(1), newR(2));
  }

  // Converte as Unidades e Limpa
  for (int i = 0; i < FockNumber; i++)
  {
    FockMatrices[i] = FockMatrices[i] * TOEV;
    FockMatrices[i].clean(UThr);
    Overlap[i]      = Overlap[i] * TOEV;
    Overlap[i].clean(UThr);
  }

  // Imprime Hamiltoniana
  cout << "# hamiltonian" << endl;
  for(int HH = 0; HH < FockNumber; HH++)
  {
    for(int i = 0; i < nOrbitals; i++)
    {  
      for(int j = 0; j < nOrbitals; j++)
      {
      cout << std::format("{: 11.8f} {: 11.8f}j    ", FockMatrices[HH](i,j), FockMatrices[HH](i,j)*0);
      }
      cout << endl; 
    } 
    cout << '&' << endl;
  }

  // Imprime Overlap                                                                               
  cout << "# overlap" << endl;                                                                      
  for(int SS = 0; SS < FockNumber; SS++)                                                                
  {                                                                                                     
    for(int i = 0; i < nOrbitals; i++)                                                                  
    {                                                                                                   
      for(int j = 0; j < nOrbitals; j++)                                                                
      {                                                                                                 
      cout << std::format("{: 11.8f} {: 11.8f}j    ", Overlap[SS](i,j), Overlap[SS](i,j)*0);  
      }                                                                                                 
      cout << endl;                                                                                     
    }                                                                                                   
    cout << '&' << endl;                                                                                
  }                                                                                                     

  // Filling
  cout << "# filling" << endl;
  cout << "4" << endl;
  cout << "#" << endl;




}

