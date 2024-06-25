#include "BandCalc.h"


int main()                                                  
{                                                           
  TightBinding hBN = TightBinding();                        
                                                            
  hBN.SetR("Cell.inp");                                     
  hBN.LoadFockMatrices(49);                                 
  hBN.LoadOverlap(49);                                      
  
  arma::mat Bands = arma::mat(8, 100, arma::fill::zeros);

  // 0.0 ---> 0.5
  double step = (0.5 - 0.0)/100;
  for(int i = 0; i < 100; i++)
    Bands.col(i) = hBN.BandCalc(step * i, 0.0, 0.0); 
                                   

  for(int i = 0; i < 8; i++)
  {
    for(int j = 0; j < 100; j++)
    {
      printf("%f\t%f\n", step * j, Bands(i, j));
    }
    printf("\n");
  }
                                                            
  return 0;                                                 
}                                                           
