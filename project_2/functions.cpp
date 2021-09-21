#include "functions.hpp"


arma::mat tridiagonal_matrix(double a, double d, double N)
{
    /*
    Fills a NxN tridiagonal matrix with elements a, b and c. (Solves problem 3)
    Args:
        tri_matrix (armadillo Matrix): matrix to be filled with values
        a          (double)          : 
        b          (double)          : 
        c          (double)          : 
        N          (int)             : 
    */

    arma::mat tri_matrix(N, N);

    for(int i=1; i<N-1; i++)
    {
        tri_matrix(i,i)  = d; // Diagonal
        tri_matrix(i,i+1) = a; // super diagonal
        tri_matrix(i,i-1) = a; // sub diagonal
    }
    // Top boundary
    tri_matrix(0,0) = d;
    tri_matrix(0,1) = a;

    // Bottom boundary
    tri_matrix(N-1,N-1) = d;
    tri_matrix(N-1,N-2) = a;

    return tri_matrix;
}

void print_array(double array[], int n)
{
    /* 
    Prints all elements in an array nicely
    Args:
        array (double array): array to be printed
        n     (int)         : length of array
    */
    for(int i = 0; i < n; i++)
    {
        std::cout << array[i] << '\n';
    }
}

void write_to_file(std::string filename, double array[], double lin[], int n)
{
    /*
    Writes the elements of two arrays, array and lin with length n, to a file with filename 'filename'
    in scientific notation
    Args:
        filename    (string) : filename of output file
        array (double array) : array whose elements will fill the second column within the file
        lin   (double array) : array whose elements will fill the first column within the file
        n     (int)          : length of array and lin
    */
    std::fstream file;
    file.open(filename+".txt", std::ios::out);

    std::cout.precision(5);

    for (int i=0; i < n; i++)
    {
        file << lin[i] << ' ' << array[i] << std::scientific << '\n';
    }
    file.close();
}

arma::vec eigenvalues(double a, double d, int n)
{
    /*
    finds eigenvalues analytically for a tridiagonal matrix with elements defined by
    a and d.
    Args:
        a (double) : 
    */
    arma::vec values(n); 

    for (int i=1; i<n; i++)
    {
        values(i) = d + 2*a*std::cos(i*3.1415/(n+1));
    }
    return values;
}
