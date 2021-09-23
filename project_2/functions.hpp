
#ifndef Functions_hpp
#define Functions_hpp

#include <armadillo>
#include <iostream>
#include <cmath>

arma::mat tridiagonal_matrix(double, double, double);
arma::vec eigenvalues(double, double, int);
arma::mat eigenvectors(int);
double largest_off_element(arma::mat, int&, int&, int);
void jacobi_rotate(arma::mat&, arma::mat&, int, int, int);
void jacobi_eigensolver(arma::mat&, arma::vec&, arma::mat&, 
                        const int, int&, bool&, int);
void linspace(arma::vec&, double, double, int);
void write_eig_to_file(arma::vec, arma::mat, int n, std::string);
void write_lin_to_file(arma::vec, int, std::string);

#endif