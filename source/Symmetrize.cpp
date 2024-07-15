#include "TBModel.h"

void TightBinding::Symmetrize()
{
  cout << ">> Begining Symmetrization Process\n";
  arma::vec Rn = arma::vec(3, arma::fill::zeros);
  arma::mat ARn;
  
  SuperCell(2);

  Rn(0) = 1;
  Rn(1) = 0;
  Rn(2) = 0;

  
  ARn = BuildNeighborn(Rn);
  cout << ARn.clean(1E-1);
  cout << "==================================\n";

  ARn = ARn*RotC3().t();
  cout << ARn.clean(1E-1);
  cout << "==================================\n";
  ARn = ARn*RotC3().t();
  cout << ARn.clean(1E-1);
}

arma::mat TightBinding::BuildNeighborn(arma::vec Rn)
{
  arma::mat ARn = arma::mat(nAtoms, 3, arma::fill::zeros);

  // Shift Cell Vector
  Rn = R.t() * Rn;

  // Shift Atoms
  for(int i = 0; i < nAtoms; i++)
  {
    ARn.row(i) = rAtoms.row(i) + Rn.t();
  }
  
  return ARn;

}

arma::vec TightBinding::CatchCell(arma::vec aRn)
{
  // shifted Coordinates
  arma::mat sRn = arma::inv(R.t());
  arma::vec Rn; 
  // do [(R^t)^(-1)] * [R^t * Rn], where Rn = (l,m,n)
  Rn = sRn * aRn;

  return arma::round(Rn);
}

arma::mat TightBinding::RotC3()
{
  double Angle = 2*PI/3;
  arma::mat RotC3 = arma::mat(3, 3);
  RotC3 = {{cos(Angle), -sin(Angle), 0},
           {sin(Angle),  cos(Angle), 0},
           {         0,           0, 1}};

  return RotC3;
}

void TightBinding::SuperCell(int N)
{
  ofstream oFile;
  oFile.open("SuperCell.xyz");

  oFile << nAtoms*(2*N+1)*(2*N+1) << endl << endl;
  arma::mat ARn;
  for(int i = -N; i <= N; i++)
  {
    for(int j = -N; j <= N; j++)
    {
      ARn = BuildNeighborn({i, j, 0});
      oFile << "B\t" << ARn.row(0);
      oFile << "N\t" << ARn.row(1);
    }
  }

  oFile.close();
}
