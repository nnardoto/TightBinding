#include "TBModel.h"


int main(int argc, char **argv)                                                  
{                        
  if(argc == 1)
  {
    cout << "Input missed. Call like \"TightBinding.x Input.dat\"" << endl;
  }
 
  // Cria o Objeto TightBinding
  TightBinding Model;         

  // Carrega a partir de um Arquivo de Input
  Model.Load(argv[1]);                                     
  Model.SymmetrizeC3();
  Model.PathCalc();
  Model.Export2Xatu();

  return 0;                                                 
}                                                           
