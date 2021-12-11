#define ARMA_DONT_USE_WRAPPER

#include <armadillo>


arma::mat calc_asymptotic_variance(arma::mat Sigma_ll, 
                                   arma::mat Sigma_ls, 
                                   arma::mat Sigma_ss,
                                   arma::field <arma::mat > Sigma_ss_blockwise,
                                   double sigma2_s, 
                                   unsigned int n,
                                   arma::mat Xl_test, 
                                   arma::mat Xs_test);

arma::mat calc_var_betal(arma::mat Sigma_ll, 
                         arma::mat Sigma_ls, 
                         arma::mat Sigma_ss,
                         arma::mat A_inverse,
                         unsigned int n);

  arma::mat calc_var_betas(arma::mat Sigma_ss, 
                           arma::mat Sigma_ls,
                           arma::mat A_inverse,
                           double sigma2_s,
                           unsigned int n,
                           arma::mat var_bl);

  arma::mat calc_A_inverse(arma::field <arma::mat > field, 
                           double sigma2_s, 
                           unsigned int n);


arma::mat BlockDiag( arma::field<arma::mat> x );

arma::mat ConcatenateColumns( arma::field<arma::mat> x );

arma::field <arma::mat> assembleMatrices(arma::field < arma::mat> field);