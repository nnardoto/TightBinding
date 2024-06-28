#include <iostream>
#include <armadillo>
#include <any>
#include <map>

#define PI 3.14159265359


using namespace std;

class TightBinding
{
  public:
    arma::vec BandCalc(double k1, double k2, double k3);
    void LoadFock();
    void LoadOverlap();
    void Load(string FileName);
    map<string, string> Parser(string FileName);

  private:
    arma::vec K[3];
    arma::vec R[3];
    int Rn[9][3];
    map<string, string> InputDict;
    int **Map;
    int Index(int i, int j, int k);
    int **aIndex;
    int FockNumber = 0;
    int NeighborsCells;
    void FockCount();
    arma::mat *FockMatrices;
    arma::mat *Overlap;
};

