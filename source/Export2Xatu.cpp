#include "TBModel.h"
#include <format>
#define TOEV  27.2114079527

void TightBinding::Export2Xatu()
{
  ofstream XatuModel;
  XatuModel.open("hBN.model");

  int Dimensions = stoi(InputDict["Dimensions"]);
  // int nOrbitals  = GetOrbitals();
  arma::mat R, Atoms;
  arma::vec newR;
  R.load(InputDict["CellFile"]);
  Atoms.load(InputDict["AtomsFile"]);
  
  XatuModel << "# dimension" << endl;
  XatuModel << Dimensions << endl;

  XatuModel << "# norbitals" << endl;
  XatuModel << format("{} {}", nOrbitals, nOrbitals) << endl;
  
  XatuModel << "# bravaislattice" << endl;
  for(int i = 0; i < Dimensions; i++)
  {
    XatuModel << std::format("{: 11.8f}    {: 11.8f}    {: 11.8f}\n", R.row(i)(0), R.row(i)(1), R.row(1)(2));
  }

  XatuModel << "# motif" << endl;
  for(int i = 0; i < 2; i++)
  {
    //TODO PADRONIZAR O INPTU
    XatuModel << std::format("{: 11.8f}    {: 11.8f}    {: 11.8f} {: 5.4f}\n", Atoms.row(i)(0), Atoms.row(i)(1), Atoms.row(i)(2), float(i));
  }

  XatuModel << "# bravaisvectors" << endl;
  for(int i = 0; i < FockNumber; i++)
  {
    arma::vec A(3);
    newR = arma::vec(3, arma::fill::zeros);
    A(0) = aIndex[i][0];
    A(1) = aIndex[i][1];
    A(2) = aIndex[i][2];
    newR = R.row(0).t()*A(0) + R.row(1).t()*A(1) + R.row(2).t()*A(2);
    XatuModel << std::format("{: 11.8f}    {: 11.8f}    {: 11.8f}\n", newR(0), newR(1), newR(2));
  }

//  // Converte as Unidades e Limpa
//  for (int i = 0; i < FockNumber; i++)
//  {
//    FockMatrices[i] = FockMatrices[i] * TOEV;
//    FockMatrices[i].clean(UThr);
//    Overlap[i] = Overlap[i] * TOEV;
//    Overlap[i].clean(UThr);
//  }

  // Imprime Hamiltoniana
  XatuModel << "# hamiltonian" << endl;
  for(int HH = 0; HH < FockNumber; HH++)
  {
    for(int i = 0; i < nOrbitals; i++)
    {  
      for(int j = 0; j < nOrbitals; j++)
      {
        XatuModel << std::format("{: 9.6f} {: 9.6f}j    ", FockMatrices[HH](i,j).real(), FockMatrices[HH](i,j).imag());
      }
      XatuModel << endl; 
    } 
    XatuModel << '&' << endl;
  }

  // Imprime Overlap                                                                               
  XatuModel << "# overlap" << endl;                                                                      
  for(int SS = 0; SS < FockNumber; SS++)                                                                
  {                                                                                                     
    for(int i = 0; i < nOrbitals; i++)                                                                  
    {                                                                                                   
      for(int j = 0; j < nOrbitals; j++)                                                                
      {                                                                                                 
        XatuModel << std::format("{: 11.7f} {: 11.7f}j    ", Overlap[SS](i,j).real(), Overlap[SS](i,j).imag());  
      }                                                                                                 
      XatuModel << endl;                                                                                     
    }                                                                                                   
    XatuModel << '&' << endl;                                                                                
  }                                                                                                     

  // Filling
  XatuModel << "# filling" << endl;
  XatuModel << "4" << endl;
  XatuModel << "#" << endl;

  // Close File
  XatuModel.close();
}

