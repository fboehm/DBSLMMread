/*
 Deterministic Bayesian Sparse Linear Mixed Model (DBSLMM)
 Copyright (C) 2019  Sheng Yang and Xiang Zhou
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <armadillo>
#include "../include/calc_asymptotic_variance.h"


using namespace std;
using namespace arma;

int main()
{
  //int nchr = 22; //number of chromosomes
  int nchr = 1;
  //initialize a arma::field to store outputs for var calcs!
  
  arma::field < arma::mat> training(nchr, 5);
  arma::field < arma::mat> test(nchr, 5);
  double sigma2_s = cPar.h / (double)cPar.nsnp;
  
  for (int i = 0; i < nchr - 1; ++i){
    int chr = i + 1;
    std::string filetr ("Chr" + std::to_string(chr) + "_training.dat");
    std::string filete ("Chr" + std::to_string(chr) + "_test.dat");
    arma::field <arma::mat > tr;
    arma::field <arma::mat > te; 
    tr.load(filetr);
    te.load(filete);
    training.row(i) = assembleMatrices(tr);//store in a two-dimensional field
    test.row(i) = assembleMatrices(te);
  }
  //var calcs here! 
  //1. assemble genome-wide matrices from "results"
  // results is a 22-long field where each entry is itself a 1d field containing 5 matrices
  arma::field < arma::mat > mats_training = assembleMatrices(training);
  arma::field < arma::mat > mats_test = assembleMatrices(test);
  //2. input matrices to calc_asymptotic_variance
  arma::mat vv = calc_asymptotic_variance(mats_training(2), 
                                          arma::trans(mats_training(1)), 
                                          mats_training(0), 
                                          cPar.h, 
                                          cPar.n, 
                                          mats_test(4), 
                                          mats_test(3));
  //3. write diagonal of var to a csv file
  arma::vec vd = diagvec(vv);
  vd.save("out.csv", csv_ascii);
  
  
  return EXIT_SUCCESS;
}
