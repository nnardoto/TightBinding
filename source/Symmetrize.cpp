#include "TBModel.h"

void TightBinding::Symmetrize()
{
  cout << ">> Begining Symmetrization Process\n";
  arma::vec Rn  = arma::vec(3, arma::fill::zeros);
  arma::vec dRn = arma::vec(3, arma::fill::zeros);
  arma::vec R0  = arma::vec(3, arma::fill::zeros);
  arma::mat ARn;
  arma::cx_mat Buffer = arma::cx_mat(nOrbitals, nOrbitals, arma::fill::zeros);
  arma::mat C3  = RotCn(3).t();
 
 int nOrb[2] = {4, 4};
 int iOrb[2] = {0, 4};


  SuperCell(3);

  int l, m, n;
  for(int ii = -3; ii <= 3; ii++)
  {
    for(int jj = -3; jj <= 3; jj++)
    {
      R0(0) = ii;
      R0(1) = jj;
      R0(2) =  0; 
      if(Index(ii, jj, 0) > 0)
      {
        for(int i = 0; i < nAtoms; i++)
        {
          Rn = R0;
          for(int j = 0; j < 3; j++)
          {
            ARn = BuildNeighborn(Rn);
            l = Rn(0); m = Rn(1); n = Rn(2);
            for(int k = 0; k < nAtoms; k++)
            {
              cout << format("{: d} {: d} --> {: d} {: d}     Line:{}-{}      Colunm:{}-{}\n", ii, jj, l, m, iOrb[k], iOrb[k] + nOrb[k] - 1, iOrb[i], iOrb[i] + nOrb[i] - 1);
             // Buffer.submat(iOrb[i], iOrb[k], iOrb[i] + nOrb[i] - 1, iOrb[k] + nOrb[k] - 1) = 
             // Buffer.submat(iOrb[i], iOrb[k], iOrb[i] + nOrb[i] - 1, iOrb[k] + nOrb[k] - 1) + 
             //  FockMatrices[Index(R0(0), R0(1), R0(2))].submat(iOrb[i], iOrb[k], iOrb[i] + nOrb[i] - 1, iOrb[k] + nOrb[k] - 1);
            }
            ARn = ARn*C3;
            ARn.clean(1E-1);
            Rn  = CatchCell(ARn.row(1).t() - rAtoms.row(1).t());
         
          }
        }
      }
    }
  }
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

arma::mat TightBinding::RotCn(int N)
{
  double Angle = 2*PI/N;
  arma::mat RotCn = arma::mat(3, 3);
  RotCn = {{cos(Angle), -sin(Angle), 0},
           {sin(Angle),  cos(Angle), 0},
           {         0,           0, 1}};

  return RotCn;
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
      ARn = BuildNeighborn({(double) i, (double)j, (double)0});
      oFile << "B\t" << ARn.row(0);
      oFile << "N\t" << ARn.row(1);
    }
  }

  oFile.close();
}
