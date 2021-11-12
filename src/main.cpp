#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>

#include "../include/dbslmm.hpp"

using namespace std;

int main(int argc, char * argv[])
{
  DBSLMM cDB;
  PARAM cPar;
  
  if (argc <= 1) {
    cDB.printHeader();
    return EXIT_SUCCESS;
  }
  if (argc == 2 && argv[1][0] == '-' && argv[1][1] == 'h') {
    cDB.printHelp();
    return EXIT_SUCCESS;
  }
  for (int i = 0; i < argc; ++i)
    cout << argv[i] << "\n";
  
  cDB.Assign(argc, argv, cPar);
  
  cDB.BatchRun(cPar);
  return EXIT_SUCCESS;
}
