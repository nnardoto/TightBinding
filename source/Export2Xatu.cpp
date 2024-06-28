#include "TBModel.h"

void TightBinding::Export2Xatu()
{
  // int Dimensions = stoi(InputData["Dimensions"]);
  // int nOrbitals  = GetOrbitals();
  // arma::vec R;
  //R.load(InputData["CellFile"]);

  // Imprime Hamiltoniana
  for(int HH = 0; HH < FockNumber; HH++)
  {
    for(int i = 0; i < 8; i++)
    {  
      for(int j = 0; j < 8; j++)
      {
        cout << FockMatrices[HH](i,j) << "  j0.000000";
      }
      cout << endl; 
    } 
    cout << '&' << endl;
  }
}

