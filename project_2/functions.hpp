
#ifndef Functions_hpp
#define Functions_hpp

#include <armadillo>
#include <iostream>
#include <cmath>

arma::mat tridiagonal_matrix(double, double, double);
arma::vec eigenvalues(double, double, int);
arma::mat eigenvectors(int);
double largest_off_element(arma::mat, int&, int&, int);

#endif