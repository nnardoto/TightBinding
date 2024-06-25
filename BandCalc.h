#include <iostream>
#include <armadillo>
#include <cstdio>

#define PI 3.14159265359


using namespace std;

class TightBinding
{
  public:
    arma::vec BandCalc(double k1, double k2, double k3);
    int LoadFockMatrices(int N);
    int LoadOverlap(int N);
    void SetR(string FileName);

  private:
    arma::vec K[3];
    arma::vec R[3];
    int Rn[9][3];
    int FockNumber = 9;
    arma::mat *FockMatrices;
    arma::mat *Overlap;
};

