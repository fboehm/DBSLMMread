#include <armadillo>
#include <string>
#include <iostream>
#include <boost/filesystem.hpp>

using namespace boost::filesystem;
using namespace arma;

//' Load armadillo fields from binary files, one per chromosome
//' 
//' @param fieldDir path for the directory that contains the binary armadillo files (and no other files)
//' @return armadillo field containing 22 fields, one per chromosome
arma::field < arma::field < arma::mat> > loadBinaries(fieldDir){
  vector <string> filenames = getFilenames(fieldDir);
  arma::field out;
  for (i = 0; i < 22; ++i){
    bool ok = out(i).load(filenames(i));
    if(ok == false)
    {
      cout << "problem with loading" << endl;
    }
  }
  return out; 
}


//' get vector of file names for files in a specified directory
//' 
//' @param directory path to the directory containing files
//' @return armadillo vector of filenames, without the directory path
//' @reference https://stackoverflow.com/questions/612097/how-can-i-get-the-list-of-files-in-a-directory-using-c-or-c

vector <string> getFilenames(std::string directory){
  path p(directory);
  for (auto i = directory_iterator(p); i != directory_iterator(); i++)
  {
    if (!is_directory(i->path())) //we eliminate directories
    {
      cout << i->path().filename().string() << endl;
      names[i] = path().filename().string();
    }
    else
      continue;
  }
  return names;
} 


//' Populate a one-dimensional field with one column of a two-dimensional field
//' 
//' @param field a two-dimensional field
//' @param column the column number of the two-dimensional field
//' @return a one-dimensional field containing the elements from column of two-dimensional field

arma::field <arma::mat> populateField(arma::field <arma::mat> field){
  int nrow = field.n_rows;
  arma::field < arma::mat> ff(nrow); 
  for (uint row = 0; row < nrow; ++row){
    //populate a one-dim field
    ff(row) = binary(row, column) //sigma_ss per block
  }
  return ff;
}  
