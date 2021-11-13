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
   
  // recall the structure of the elements of binaries. Each entry is itself a field, and those 
  // fields each have five entries. 
  arma::field <arma::mat> Sigma_ss(22);
  arma::field <arma::mat> Sigma_sl(22);
  arma::field <arma::mat> Sigma_ll(22);
  arma::field <arma::mat> geno_s(22);
  arma::field <arma::mat> geno_l(22);
  
  //assemble chromosome-specific matrices
  for (uint i = 0; i < 22; ++i){
    string chr = to_string(i + 1);
    string filename = "binary" + chr  + ".dat";
    arma::field < arma::mat> binary = loadBinary(filename);
    arma::field <arma::mat> onecol = populateField(binary, 0);
    arma::mat Sigma_ss(i) = BlockDiag(onecol);
    arma::field <arma::mat> onecol = populateField(binary, 1);
    arma::mat Sigma_sl(i) = BlockDiag(onecol);
    arma::field <arma::mat> onecol = populateField(binary, 2);
    arma::mat Sigma_ll(i) = BlockDiag(onecol);
    arma::field <arma::mat> onecol = populateField(binary, 3);
    arma::mat geno_ss(i) = ConcatenateColumns(onecol);
    arma::field <arma::mat> onecol = populateField(binary, 4);
    arma::mat geno_ll(i) = ConcatenateColumns(onecol);
  }
  // combine 5 matrices into genome-wide versions
  arma::mat geno_ss_matrix = ConcatenateColumns(geno_ss);
  arma::mat geno_ll_matrix = ConcatenateColumns(geno_ll);
  arma::mat Sigma_ll_matrix = BlockDiag(Sigma_ll);
  arma::mat Sigma_sl_matrix = BlockDiag(Sigma_sl);
  arma::mat Sigma_ss_matrix = BlockDiag(Sigma_ss);
  // calculate asympt var
  arma::mat vv = calc_asymptotic_variance(Sigma_ll_matrix, 
                                          Sigma_sl_matrix.t(), 
                                          Sigma_ss_matrix, 
                                          sigma2_s,
                                          n, 
                                          Xl_test, Xs_test);
  vv.arma::save("var.csv", csv_ascii);
  return EXIT_SUCCESS;
}
