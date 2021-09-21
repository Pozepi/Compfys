#include "functions.hpp"

int main()
{
    double N = 6;
    double h = 1/(N+1);
    double a = -1/(h*h);
    double d =  2/(h*h);

    arma::mat tri_matrix = tridiagonal_matrix(a, d, N);
    std::cout << tri_matrix << '\n';

    //arma::vec eigval = arma::eig_sym(tri_matrix);
    /*arma::mat eigvec;
    
    arma::eig_sym(eigval, eigvec, tri_matrix);
    */
    //std::cout << eigval << '\n';
    //std::cout << eigvec << '\n';
    

    arma::vec eigenval = eigenvalues(a, d, N);
    std::cout << eigenval;
    
    return 0;
}





