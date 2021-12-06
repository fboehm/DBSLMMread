#include <armadillo>
#include <string>
#include <iostream>
#include <boost/filesystem.hpp>


arma::field < arma::field < arma::mat> > loadBinaries(fieldDir);

vector <string> getFilenames(std::string directory);

arma::field <arma::mat> populateField(arma::field <arma::mat> field);

