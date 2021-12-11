#define ARMA_DONT_USE_WRAPPER

#include <armadillo>


arma::mat calc_asymptotic_variance(const arma::mat& Sigma_ll, 
                                   const arma::mat& Sigma_ls, 
                                   const arma::mat& Sigma_ss,
                                   const arma::field <arma::mat >& Sigma_ss_blockwise,
                                   double sigma2_s, 
                                   unsigned int n,
                                   const arma::mat& Xl_test, 
                                   const arma::mat& Xs_test);

arma::mat calc_var_betal(const arma::mat& Sigma_ll, 
                         const arma::mat& Sigma_ls, 
                         const arma::mat& Sigma_ss,
                         const arma::mat& A_inverse,
                         unsigned int n);

  arma::mat calc_var_betas(const arma::mat& Sigma_ss, 
                           const arma::mat& Sigma_ls,
                           const arma::mat& A_inverse,
                           double sigma2_s,
                           unsigned int n,
                           const arma::mat& var_bl);

  arma::mat calc_A_inverse(const arma::field <arma::mat >& field, 
                           double sigma2_s, 
                           unsigned int n);


arma::mat BlockDiag( const arma::field<arma::mat>& x );

arma::mat ConcatenateColumns( const arma::field<arma::mat>& x );

arma::field <arma::mat> assembleMatrices(const arma::field < arma::mat>& field);