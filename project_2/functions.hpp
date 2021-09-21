
#ifndef Functions_hpp
#define Functions_hpp

#include <armadillo>
#include <iostream>
#include <cmath>

arma::mat tridiagonal_matrix(double, double, double);
arma::vec eigenvalues(double, double, int);

#endif