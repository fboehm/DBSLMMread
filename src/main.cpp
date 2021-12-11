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
#include "../include/calc_asymptotic_variance.hpp"
#include "../include/dbslmm.hpp"

using namespace std;
using namespace arma;

int main(int argc, char * argv[])
{
  DBSLMM cDB;
  PARAM cPar;
  cDB.Assign(argc, argv, cPar);
  
  //int nchr = 22; //number of chromosomes
  int nchr = 1;
  //initialize a arma::field to store outputs for var calcs!
  
  arma::field < arma::mat> training(5000, 5); //5000 is always bigger than the number of LD BLocks in the genome
  arma::field < arma::mat> test(5000, 5);
  double sigma2_s = cPar.h / (double)cPar.nsnp;
  unsigned int row_total = 0; //initialize counter for number of rows in field
  for (int i = 0; i < nchr; ++i){
    //int i = 0;
    int chr = i + 1;
    std::string filetr ("Chr" + std::to_string(chr) + "_training.dat");
    std::string filete ("Chr" + std::to_string(chr) + "_test.dat");
    arma::field <arma::mat > tr;
    arma::field <arma::mat > te; 
    tr.load(filetr);
    cout << filetr << endl;
    cout << "number of elements in tr: " << tr.n_elem << endl;
    cout << "number of rows in tr: " << tr.n_rows << endl;
    cout << "number of columns in tr: " << tr.n_cols << endl;
    
    te.load(filete);
    cout << filete << endl; 
    cout << "number of elements in te: " << te.n_elem << endl;
    cout << "number of rows in te: " << te.n_rows << endl;
    cout << "number of columns in te: " << te.n_cols << endl;
    
    unsigned int rows_in_block = te.n_rows;
    training.rows(row_total, row_total+rows_in_block - 1) = tr;
    test.rows(row_total, row_total+rows_in_block - 1) = te;
    row_total = row_total + rows_in_block;
    
//    training.row(i) = assembleMatrices(tr);//store in a two-dimensional field
//    test.row(i) = assembleMatrices(te);
  }
  training = training.rows(0, row_total - 1);
  cout << "training has this many rows: " << row_total << endl;
  test = test.rows(0, row_total - 1);
  cout << "test has this many rows: " << row_total << endl;  
  //var calcs here! 
  //1. assemble genome-wide matrices from "training" & "test"
  arma::field < arma::mat > mats_training = assembleMatrices(training);
  arma::field < arma::mat > mats_test = assembleMatrices(test);
  //2. input matrices to calc_asymptotic_variance
  cout << "Starting asymptotic variance calculations..." << endl; 
  arma::mat vv = calc_asymptotic_variance(mats_training(2), //Sigma_ll
                                          arma::trans(mats_training(1)), // Sigma_ls 
                                          mats_training(0), //Sigma_ss
                                          training.col(0),
                                          sigma2_s, 
                                          cPar.n, 
                                          mats_test(4), //X_l
                                          mats_test(3)); // X_s
  //3. write diagonal of var to a csv file
  arma::vec vd = diagvec(vv);
  vd.save("out.csv", csv_ascii);

  return EXIT_SUCCESS;
}
