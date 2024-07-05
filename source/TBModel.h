#include <iostream>
#include <armadillo>
#include <any>
#include <map>

#define PI 3.14159265359


using namespace std;

class TightBinding
{
  public:
    void PathCalc();
    void Export2Xatu();
    void Load(string FileName);
    map<string, string> Parser(string FileName);
    arma::cx_mat GetH(double k1, double k2, double k3);

  private:
    arma::vec K[3];
    arma::vec R[3];

    // Used for indexing matrieces and vectors
    int **Map;
    int **aIndex;
    double UThr = 1.0E-13;   


    // Used for iterations
    int FockNumber       = 0;
    int NeighborsCells   = 0;
    int nOrbitals        = 0;
    int WanN             = 0;
    bool OrthogonalBasis = true;

    // Store input Information
    map<string, string> InputDict;
    void MakeMap(int Neighbors);

    // Matrices e Degenerecencia
    arma::cx_mat *FockMatrices;
    arma::cx_mat *Overlap;
    arma::vec *Degen;

    // Private Mathods
    int Index(int i, int j, int k);
    void FockCount();
    void LoadFock();   
    void LoadOverlap();
    void WF_SkipHead(ifstream &iFile);
    void OMX_Load();
    void W90_Load();

    void WF_LoadFock();
    arma::vec    BandCalc(double k1, double k2, double k3);
    arma::vec WF_BandCalc(double k1, double k2, double k3);
};

