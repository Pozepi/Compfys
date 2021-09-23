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

    for (int i=1; i<n+1; i++)
    {
        values(i-1) = d + 2*a*std::cos(i*3.1415/(n+1));
    }
    return values;
}

arma::mat eigenvectors(int n)
{
    /* DOCSTRING */
    arma::mat vectors(n,n);
    
    for (int i=1; i<n+1; i++)
    {
        for (int j=1; j<n+1; j++)
        {
            vectors(i-1,j-1) = std::sin(j*i*3.1415/(n+1));
        }
    }
    arma::mat transposed = vectors.t();
    return transposed;
}

double largest_off_element(arma::mat A, int& k, int& l, int n)
{
    /*
    Finds largest element value and index off of the diagonal in a symetric matrix
    Args:
        A (armada matrix) : Matrix to search for the max element in 
        k (int&)          : used to record index of max element
        l (int&)          : used to record index of max element
        n (int)           : length of the matrix (matrix is nxn)
    */
    double max = -1;
    int null;
    for (int i=0; i<n; i++)
    {
        for (int j=i+1; j<n; j++)
        {
            if (i != j)
            {
                if (fabs(A(i,j)) > max)
                {   
                    max = fabs(A(i,j));
                    l = i;
                    k = j;
                }
            }
        }
    }
    return max;
}

// Performs a single Jacobi rotation, to "rotate away"
// the off-diagonal element at A(k,l).
// - Assumes symmetric matrix, so we only consider k < l
// - Modifies the input matrices A and R
void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l, int n)
{
    //double theta = std::arctan(phi + sqrt(phi*phi + 1));
    //for (int i=0; i>n;i++)
    
    double c;
    double s;
    double phi = (A(l,l) - A(k,k))/(2*A(k,l));

    double t = 1/(fabs(phi) + sqrt(phi*phi+1)); 

    if(phi<0)
        t = -t;  
    
    c = 1/(sqrt(1+t*t)); // cosine
    s = c*t;             // sine

    //arma::mat tmp = A; // bii = aii included

    // b_kk = a_kk * cos^2θ − 2 * a_kl * cosθ * sinθ + a_ll * sin^2θ
    /*
    A(k,k) = tmp(k,k)*c*c - 2.0*tmp(k,l)*c*s + tmp(l,l)*s*s;
    A(l,l) = tmp(k,k)*s*s + 2.0*tmp(k,l)*c*s + tmp(l,l)*c*c;
    A(k,l) = (A(k,k) - A(l,l))*c*s + A(k,l)*(c*c - s*s);
    //A(l,k) = 0.0;

    for(int i = 0; i<n; i++)
    {
        if (i != k && i != l)
        {
            // b_ik = a_ik * cosθ − a_il * sinθ, i != k, i != l 
            // b_il = a_il * cosθ + a_ik * sinθ, i != k, i != l
            A(i,k) = tmp(i,k)*c - tmp(i,l)*s;
            A(i,l) = tmp(i,l)*c + tmp(i,k)*s;
            A(k,i) = tmp(i,k);
            A(l,i) = tmp(i,l);
        }

        double tmp_ik = R(i,k);
        double tmp_il = R(i,l);
        R(i,k) = tmp_ik*c - tmp_il*s;
        R(i,l) = tmp_il*c + tmp_ik*s;
    }
    */
    A(k,k) = A(k,k)*c*c - 2.0*A(k,l)*c*s + A(l,l)*s*s;
    A(l,l) = A(k,k)*s*s + 2.0*A(k,l)*c*s + A(l,l)*c*c;
    A(k,l) = 0.0; //(A(k,k) - A(l,l))*c*s + A(k,l)*(c*c - s*s);
    A(l,k) = 0.0;

    for(int i = 0; i<n; i++)
    {
        if (i != k && i != l)
        {
            // b_ik = a_ik * cosθ − a_il * sinθ, i != k, i != l 
            // b_il = a_il * cosθ + a_ik * sinθ, i != k, i != l
            A(i,k) = A(i,k)*c - A(i,l)*s;
            A(i,l) = A(i,l)*c + A(i,k)*s;
            A(k,i) = A(i,k);
            A(l,i) = A(i,l);
        }

        double tmp_ik = R(i,k);
        double tmp_il = R(i,l);
        R(i,k) = tmp_ik*c - tmp_il*s;
        R(i,l) = tmp_il*c + tmp_ik*s;
    }
}

// Jacobi method eigensolver:
// - Runs jacobo_rotate until max off-diagonal element < eps
// - Writes the eigenvalues as entries in the vector "eigenvalues"
// - Writes the eigenvectors as columns in the matrix "eigenvectors"
//   (The returned eigenvalues and eigenvectors are sorted using arma::sort_index)
// - Stops if it the number of iterations reaches "maxiter"
// - Writes the number of iterations to the integer "iterations"
// - Sets the bool reference "converged" to true if convergence was reached before hitting maxiter
void jacobi_eigensolver(arma::mat& A, arma::vec& eigenvalues, arma::mat& eigenvectors, 
                        const int maxiter, int& counter, bool& converged, int n)
{
    double eps = 1e-15;
    double sum = 0;
    int k;
    int l; 
    while (!converged)
    {
        counter++;
        double max = largest_off_element(A, k, l, n);
        // std::cout << "k : " << k << " l : " << l << "\n";
        // std::cout << "iteration :" << counter << "\n" << A << "\n";
        jacobi_rotate(A, eigenvectors, k, l, n);

        double sum = 0;
        for(int i = 0; i<n-1; i++)
        {
            for(int j = i+1; j<n; j++)
            {
                sum += fabs(A(i,j));
            }
        }
        if (sum < eps || counter > maxiter)
        {
            converged = true;
        }
    } 
    if (sum < eps)
    {
        for (int i=0;i<n;i++)
        {
            eigenvalues(i) = A(i,i);
        }
        std::cout << "SUCEED \n";
    }
    if (counter >= maxiter && sum > eps)
    {
        std::cout << "Failed \n";
    }
}