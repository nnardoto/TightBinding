#include <iostream>
#include <armadillo>
#include <any>
#include <map>

#define PI 3.14159265359


using namespace std;

class TightBinding
{
  public:
    arma::mat PathCalc(arma::vec Origin, arma::vec Destination, int N);
    void Export2Xatu();
    void Load(string FileName);
    map<string, string> Parser(string FileName);

  private:
    arma::vec K[3];
    arma::vec R[3];

    // Used for indexing matrieces and vectors
    int **Map;
    int **aIndex;
    


    // Used for iterations
    int FockNumber = 0;
    int NeighborsCells;
    int nOrbitals = 8;

    // Store input Information
    map<string, string> InputDict;

    // Matrices 
    arma::mat *FockMatrices;
    arma::mat *Overlap;

    // Private Mathods
    int Index(int i, int j, int k);
    void FockCount();
    void LoadFock();    
    void LoadOverlap(); 
    arma::vec BandCalc(double k1, double k2, double k3);
};

