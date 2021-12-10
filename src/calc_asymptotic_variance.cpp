#define ARMA_DONT_USE_WRAPPER

#include <armadillo>

#include "../include/calc_asymptotic_variance.hpp"

using namespace std;
using namespace arma;

//' Calculate the asymptotic variance for the predicted y values
//' 
//' @param Sigma_ll Sigma_ll matrix for the whole genome
//' @param Sigma_ls Sigma_ls matrix for the whole genome
//' @param Sigma_ss Sigma_ss matrix for the whole genome
//' @param sigma2_s estimated value of sigma^2_s
//' @param n sample size for observed data (not the reference panel)
//' @param Xl_test genotypes matrix for large effect SNPs for test subjects
//' @param Xs_test genotypes matrix for small effect SNPs for test subjects
//' @return variance of predicted y values

arma::mat calc_asymptotic_variance(arma::mat Sigma_ll, 
                                   arma::mat Sigma_ls, 
                                   arma::mat Sigma_ss, 
                                   double sigma2_s, 
                                   unsigned int n,
                                   arma::mat Xl_test, 
                                   arma::mat Xs_test){
  arma::mat Ainv = calc_A_inverse(Sigma_ss, sigma2_s, n);
  arma::mat var_bl = calc_var_betal(Sigma_ll, 
                                    Sigma_ls, 
                                    Sigma_ss, 
                                    Ainv, 
                                    n);
  arma::mat var_bs = calc_var_betas(Sigma_ss, 
                                    Sigma_ls,
                                    Ainv,
                                    sigma2_s,
                                    n,
                                    var_bl);
  arma::mat result = Xl_test * var_bl * arma::trans(Xl_test) + Xs_test * var_bs * arma::trans(Xs_test);
  return(result);
}

//' Calculate A inverse matrix
//' 
//' @details (sigma^{-2}n^{-1} I_ms + Sigma_ss) = A.
//' @param Sigma_ss Sigma_ss matrix for a single LD block
//' @param sigma2_s estimate of sigma^2_s
//' @param n sample size
//' @return A inverse matrix

arma::mat calc_A_inverse(arma::mat Sigma_ss, 
                         double sigma2_s, 
                         unsigned int n){
  //determine m_s
  unsigned int m_s = Sigma_ss.n_rows;
  // calculate result
  arma::mat result = arma::inv_sympd(arma::eye(m_s, m_s) / (n * sigma2_s) + Sigma_ss);
  return result;
}



//' Calculate variance of coefficient estimator for large effects
//' 
//' @param Sigma_ll Sigma_ll constructed for one LD block
//' @param Sigma_ls Sigma_ls constructed for one LD block
//' @param Sigma_ss Sigma_ss constructed for one LD block
//' @param A_inverse inverse of (sigma^{-2}n^{-1} I_ms + Sigma_ss)
//' @param n sample size
//' @return covariance matrix

arma::mat calc_var_betal(arma::mat Sigma_ll, 
                      arma::mat Sigma_ls, 
                      arma::mat Sigma_ss,
                      arma::mat A_inverse,
                      unsigned int n){
  //calculate second matrix
  arma::mat big = Sigma_ll - Sigma_ls * A_inverse * arma::trans(Sigma_ls);
  //invert and divide by n
  arma::mat result = arma::inv_sympd(big) / n;
  return (result);
}

//' Calculate variance of coefficient estimator for small effects
//' 
//' @param Sigma_ss Sigma_ss matrix 
//' @param Sigma_ls Sigma_ls matrix 
//' @param A_inverse A inverse matrix 
//' @param sigma2_s estimated value of sigma^2_s
//' @param n sample size
//' @param var_bl variance of beta hat l
//' @return covariance matrix
  
arma::mat calc_var_betas(arma::mat Sigma_ss, 
                         arma::mat Sigma_ls,
                         arma::mat A_inverse,
                         double sigma2_s,
                         unsigned int n,
                         arma::mat var_bl){
  arma::mat small = arma::trans(Sigma_ls) - Sigma_ss * A_inverse * arma::trans(Sigma_ls);
  arma::mat term2 = small * var_bl * arma::trans(small);
  arma::mat term1 = Sigma_ss - Sigma_ss * A_inverse * Sigma_ss;
  arma::mat result = sigma2_s * sigma2_s * n * (term1 + term2);
  return (result);
}


//' Construct a block diagonal matrix from a collection of matrices
//' 
//' @param x a field of matrices, possibly of different sizes. Some matrices may have no rows and no columns
//' @return a block diagonal matrix
//' @reference https://stackoverflow.com/questions/29198893/block-diagonal-matrix-armadillo

arma::mat BlockDiag( arma::field<arma::mat> x ) {
  
  unsigned int len = x.n_elem;
  int drow = 0;
  int dcol = 0;
  arma::ivec rvec(len);
  arma::ivec cvec(len);
  //get dimensions of each matrix in the field
  for(unsigned int i = 0; i < len; i++) {
    rvec(i) = x(i).n_rows ; 
    cvec(i) = x(i).n_cols ; 
    drow += rvec(i);
    dcol += cvec(i);
  }
  //initialize matrix to be returned
  arma::mat X(drow, dcol, fill::zeros);
  int idx_row = 0;
  int idx_col = 0;
  // place matrices at correct places
  for(unsigned int i=0; i < len; i++) {
    if (rvec(i) > 0 && cvec(i) > 0){
      X.submat(idx_row, 
               idx_col, 
               idx_row + rvec(i) - 1, 
               idx_col + cvec(i) - 1) = x(i) ;
      idx_row = idx_row + rvec(i) ;
      idx_col = idx_col + cvec(i);
    }
  }
  return(X);
}


//' Construct a n by p matrix from field containing matrices with n rows, but possibly fewer columns.
//' 
//' @param x a field of matrices, possibly of different sizes, but all with the same number of rows. 
//' @return a matrix

arma::mat ConcatenateColumns( arma::field<arma::mat> x ) {
  
  unsigned int len = x.n_elem;
  //unsigned int nrow = x(1).n_rows;//problem!
  int dcol = 0;
  
  arma::ivec cvec(len);
  arma::ivec rvec(len);
  //get number of columns of each matrix in the field
  for(unsigned int i = 0; i < len; i++) {
    cvec(i) = x(i).n_cols ; 
    rvec(i) = x(i).n_rows;
    dcol += cvec(i);
  }
  unsigned int nrow = max(rvec);
  //initialize matrix to be returned
  arma::mat X(nrow, dcol, fill::zeros);
  cout << "number of rows: " << nrow << endl; 
  cout << "number of columns: " << dcol << endl; 
  
  int idx_col = 0;
  // place matrices at correct places
  for(unsigned int i=0; i < len; i++) {
    if (cvec(i) > 0){
      X.submat(0, 
               idx_col, 
               nrow - 1, 
               idx_col + cvec(i) - 1) = x(i) ;
      idx_col = idx_col + cvec(i);
    }
  }
  return(X);
}


//' Assemble one set of five matrices for one chromosome
//' 
//' @details Input is an two-dimensional arma::field, say from one chromosome, where each cell contains an arma::mat
//'     Specifically, it is a k by 5 arma::field, where k is the number of blocks on the chromosome of interest.
//' @param field a two-dimensional arma::field. See details.
//' @return a one-dimensional field containing exactly five arma::mat matrices: Sigma_ss, Sigma_sl, Sigma_ll, geno_s, geno_l    

arma::field <arma::mat> assembleMatrices(arma::field < arma::mat> field){
  arma::field <arma::mat> result(5);
  result(0) = BlockDiag(field.col(0));
  result(1) = BlockDiag(field.col(1));
  result(2) = BlockDiag(field.col(2));
  result(3) = ConcatenateColumns(field.col(3));
  result(4) = ConcatenateColumns(field.col(4));
  return result;
} 
