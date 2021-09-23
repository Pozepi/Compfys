#include "functions.hpp"

int main()
{
    double N = 3;
    double h = 1/(N+1);
    //double a = -1/(h*h);
    //double d =  2/(h*h);

    double a = 1;
    double d = -0.7;
    arma::mat tri_matrix = tridiagonal_matrix(a, d, N);
    //std::cout << tri_matrix << '\n';

    //arma::vec eigval = arma::eig_sym(tri_matrix);
    //arma::mat eigvec;
    
    //arma::eig_sym(eigval, eigvec, tri_matrix);
    
    //std::cout << eigval << '\n';
    //std::cout << eigvec << '\n';

    arma::vec eigenval = eigenvalues(a, d, N);
    arma::mat eigenvec = eigenvectors(N);
    arma::mat normalised = arma::normalise(eigenvec);
    //std::cout << "-----Analytical----- \n";
    //std::cout << eigenval << "\n";
    //std::cout << normalised << "\n";

    //std::cout << "-----Max Element----- \n";
    arma::mat test_matrix = {{1, 0, 0, 0.5}, {0, 1, -0.7, 0}, {0, -0.7, 1, 0}, {0.5, 0, 0, 1}};
    arma::mat test_matrix_2 = {{1, sqrt(2), 2}, {sqrt(2), 3, sqrt(2)}, {2, sqrt(2), 1}};
    int k;
    int l;
    double max_element = largest_off_element(test_matrix, k, l, 4);
    //std::cout << "Largest off element: " << max_element << "\n";
    //std::cout << "Index of largest off element: " << k << " " << l << "\n";

    //std::cout << "-----Rotation----- \n";
    int counter = 0;

    //std::cout << test_matrix << '\n';
    arma::mat R;
    R.eye(N,N);

    std::cout << R << "\n";
    arma::vec V(N);
    int maxiter = N*N*N;
    bool converged = false;
    std::cout << test_matrix<< "\n";
    jacobi_eigensolver(test_matrix_2, V, R, maxiter, counter, converged, N);
    std::cout << test_matrix_2<< "\n";
    std::cout << R << "\n";
    std::cout << V << "\n";
    std::cout << counter << "\n";

    /*
    std::cout << "------ARMADILLO------ \n";
    test_matrix = {{1, 0, 0, 0.5}, {0, 1, -0.7, 0}, {0, -0.7, 1, 0}, {0.5, 0, 0, 1}};
    arma::vec eigval = arma::eig_sym(test_matrix);
    arma::mat eigvec;
    
    arma::eig_sym(eigval, eigvec, test_matrix);
    
    std::cout << eigval << '\n';
    std::cout << eigvec << '\n';
    */
    return 0;
}


//g++ main.cpp functions.cpp -o test.out -O2 -larmadillo -llapack -lblas



