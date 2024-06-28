#include "TBModel.h"
#include <map>
#include <string>

map<string, string> TightBinding::Parser(string FileName)
{
  bool toGet  = true;
 
  // Dados para salvar input em forma de dicionário 
  map<string, string> InputData;

  // Arquivo de Input
  ifstream InputFile(FileName);

  // Lê arquivo linha a linha
  for(std::string Line; getline(InputFile, Line);)
  {
    stringstream ss(Line);
    string key, element;

    ss >> key;
    ss >> element;
    InputData.insert({key, element});
  }

  // Fecha arquivo
  InputFile.close();

  return InputData;
}
