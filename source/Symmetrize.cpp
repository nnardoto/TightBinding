#include "TBModel.h"

void TightBinding::Symmetrize()
{
  cout << ">> Begining Symmetrization Process\n";
  arma::vec Rn  = arma::vec(3, arma::fill::zeros);
  arma::vec R0  = arma::vec(3, arma::fill::zeros);
  arma::vec Rc  = arma::vec(3, arma::fill::zeros);
  arma::mat ARn;
  arma::mat C3  = RotCn(3).t();
  Rc = R.row(1).t()/2.0;
 
 int nOrb[2] = {4, 4};
 int iOrb[2] = {0, 4};


  SuperCell(3);

  int l, m, n;
  for(int ii = -1; ii <= 1; ii++)
  {
    for(int jj = -1; jj <= 1; jj++)
    {
      R0(0) = ii;
      R0(1) = jj;
      R0(2) =  0; 
      if(Index(ii, jj, 0) > 0)
      {
        for(int i = 0; i < nAtoms; i++)  // Atomo de Referencia, tenho verificar com calma
        {
          for(int j = 0; j < 3; j++) // C3
          Rn = R0;
          {
            ARn = BuildNeighborn(Rn);

            l = Rn(0); m = Rn(1); n = Rn(2);
            if(Index(l, m, n) > 0)
            {
              for(int k = 0; k < nAtoms; k++)
              {
                cout << "---------------------------------------------------\n";
                cout << format("{: d} {: d} 0 --> {: d} {: d} 0     Line:{}-{}      Colunm:{}-{}\n", ii, jj, l, m, iOrb[i], iOrb[i] + nOrb[i] - 1, iOrb[k], iOrb[k] + nOrb[k] - 1);
                cout << "---------------------------------------------------\n";
                cout << real(FockMatrices[Index(ii, jj, 0)].submat(iOrb[i], iOrb[k], iOrb[i] + nOrb[i] - 1, iOrb[k] + nOrb[k] - 1)); 
                cout << "---------------------------------------------------\n";
                cout << real(FockMatrices[Index( l,  m, 0)].submat(iOrb[i], iOrb[k], iOrb[i] + nOrb[i] - 1, iOrb[k] + nOrb[k] - 1)); 
                cout << "---------------------------------------------------\n";
             // Buffer.submat(iOrb[i], iOrb[k], iOrb[i] + nOrb[i] - 1, iOrb[k] + nOrb[k] - 1) + 
             //  FockMatrices[Index(R0(0), R0(1), R0(2))].submat(iOrb[i], iOrb[k], iOrb[i] + nOrb[i] - 1, iOrb[k] + nOrb[k] - 1);
              }
            }
            cout << "###################################################\n";
            ARn = ARn*C3;
            ARn.clean(1E-1);
            Rn  = CatchCell(ARn.row(0).t() - rAtoms.row(0).t() + Rc);
          }
        }
      }
      cout << "========================================\n";
    }
  }
}

void TightBinding::SymmetrizeC3()
{
  cout << ">> Begining Symmetrization Process\n";  
  arma::vec Rn  = arma::vec(3, arma::fill::zeros); 
  arma::vec R0  = arma::vec(3, arma::fill::zeros); 
  arma::vec Rc  = arma::vec(3, arma::fill::zeros); 
  arma::cx_mat buffer = arma::cx_mat(nOrbitals, nOrbitals, arma::fill::zeros);
  arma::cxmat SymFock[35];
  arma::mat ARn;                                   
  arma::mat C3  = RotCn(3).t();                    
  Rc = R.row(1).t()/2.0;                           
                                                   
  int nOrb[2] = {4, 4};                             
  int iOrb[2] = {0, 4};                             

  int ToChange[3][3];

  int n = 0;
  for(int l = -1; l <= 1; l++)
  {
    for(int m = -1; m < 1; m++)
    {
      Rn = {(double)l, (double)m, (double)n};
      ARn = BuildNeighborn(Rn);

      for(int atom = 1; atom < 2; atom++)
      {
        for(int k = 0; k < 3; k++)    // C3 Rotation
        {
          int ll = Rn(0); 
          int mm = Rn(1); 
          int nn = Rn(2); 

          cout << format("{} {} {} --> {} {} {}\n", l, m, n, ll, mm, nn);
          cout << "====================================================\n"; 
          cout << real(FockMatrices[Index( l,  m, 0)].submat(iOrb[atom], iOrb[atom], arma::size(nOrb[atom], nOrb[atom])));
          cout << "----------------------------------------------------\n"; 
          cout << real(FockMatrices[Index(ll, mm, 0)].submat(iOrb[atom], iOrb[atom], arma::size(nOrb[atom], nOrb[atom])));
          cout << "====================================================\n"; 
          
         // buffer.submat(iOrb[atom], iOrb[atom], arma::size(nOrb[atom], nOrb[atom]))
         // += FockMatrices[Index(ll, mm, 0)].submat(iOrb[atom], iOrb[atom], arma::size(nOrb[atom], nOrb[atom]))/3.0;
          
         // ToChange[k][0] = ll; ToChange[k][1] = mm; ToChange[k][2] = nn;

          ARn = ARn * C3;                                      
          ARn.clean(1E-1);                                     
          Rn  = CatchCell(ARn.row(1).t() - rAtoms.row(1).t()); 
        }

        cout << "####################################################\n";
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
