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
    double c;
    double s;
    double phi = (A(l,l) - A(k,k))/(2*A(k,l));

    double t = 1/(fabs(phi) + sqrt(phi*phi+1)); 

    if(phi<0)
        t = -t;
    
    c = 1/(sqrt(1+t*t)); // cosine
    s = c*t;             // sine

    double akk = A(k,k);
    
    A(k,k) = akk*c*c - 2.0*A(k,l)*c*s + A(l,l)*s*s;
    A(l,l) = akk*s*s + 2.0*A(k,l)*c*s + A(l,l)*c*c;
    A(k,l) = 0.0;
    A(l,k) = 0.0;

    for(int i = 0; i<n; i++)
    {
        if (i != k && i != l)
        {
            double aik = A(i,k);
            A(i,k) = aik*c - A(i,l)*s;
            A(i,l) = A(i,l)*c + aik*s;
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
    double eps = 1e-6;
    double sum = 0;
    int k;
    int l; 
    while (!converged)
    {
        counter++;
        double max = largest_off_element(A, k, l, n);
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

void write_eig_to_file(arma::vec V, arma::mat R, int n, std::string filename)
{
    std::fstream file;
    file.open(filename+".txt", std::ios::out);

    R = R.t();
    for (int i=0; i < n; i++)
    {
        file << V(i) << '\n';
        for (int j=0; j < n; j++)
        {
            file << R(i,j) << '\n';
        }
        file << '\n';
    }
    file.close();
}

void write_lin_to_file(arma::vec x, int n, std::string filename)
{
    std::fstream file;
    file.open(filename+".txt", std::ios::out);
    for (int i=0; i < n; i++)
    {
        file << x(i) << '\n';
    }
    file.close();
}
#include "functions.hpp"


arma::mat tridiagonal_matrix(double a, double d, double N)
{
    /*
    Returns a NxN tridiagonal matrix with elements d along diagonal and a adjasent (Solves problem 3)
    Args:
        a          (double)          : values adjasent to diagonal 
        b          (double)          : diagonal
        N          (int)             : size of "side" of array
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
    Finds eigenvalues analytically for a tridiagonal matrix with elements defined by
    a and d.
    Args:
        a (double) : elements adjasent to diagonal
        d (double) : diagonal
        n (int)    : size of vector
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
    /* 
    Finds the eigenvectors for matrix A = diag(a,d,a)
    Args:
        n (int) : size of side of matrix
    */
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
        k (int&)          : used to record index of max element (along y) 
        l (int&)          : used to record index of max element (along x)
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
    /*
    Performs a single Jacobi rotation, to "rotate away" the off-diagonal element at A(k,l).
    - Assumes symmetric matrix, so we only consider k < l
    - Modifies the input matrices A and R
    Args:
        A (armata matrix&) : Matrix to be treated by algorithm
        R (armada matrix&) : Matrix to fill with eigenvectors
        k (int)            : y position of max elements
        l (int)            : x position of max elements
        n (int)            : size of matrix (sides, acutally nxn)
    */
    double c;
    double s;
    double phi = (A(l,l) - A(k,k))/(2*A(k,l));

    double t = 1/(fabs(phi) + sqrt(phi*phi+1)); 

    if(phi<0)
        t = -t;
    
    c = 1/(sqrt(1+t*t)); // cosine
    s = c*t;             // sine

    double akk = A(k,k);
    
    A(k,k) = akk*c*c - 2.0*A(k,l)*c*s + A(l,l)*s*s;
    A(l,l) = akk*s*s + 2.0*A(k,l)*c*s + A(l,l)*c*c;
    A(k,l) = 0.0;
    A(l,k) = 0.0;

    for(int i = 0; i<n; i++)
    {
        if (i != k && i != l)
        {
            double aik = A(i,k);
            A(i,k) = aik*c - A(i,l)*s;
            A(i,l) = A(i,l)*c + aik*s;
            A(k,i) = A(i,k);
            A(l,i) = A(i,l);
        }

        double tmp_ik = R(i,k);
        double tmp_il = R(i,l);
        R(i,k) = tmp_ik*c - tmp_il*s;
        R(i,l) = tmp_il*c + tmp_ik*s;
    }
}

void jacobi_eigensolver(arma::mat& A, arma::vec& eigenvalues, arma::mat& eigenvectors, 
                        const int maxiter, int& counter, bool& converged, int n)
{
    /*
    Jacobi method eigensolver:
    - Runs jacobo_rotate until max off-diagonal element < eps
    - Writes the eigenvalues as entries in the vector "eigenvalues"
    - Writes the eigenvectors as columns in the matrix "eigenvectors"
    (The returned eigenvalues and eigenvectors are sorted using arma::sort_index)
    - Stops if it the number of iterations reaches "maxiter"
    - Writes the number of iterations to the integer "iterations"
    - Sets the bool reference "converged" to true if convergence was reached before hitting maxiter
    Args:
        A (armada matrix&)            : Matrix to be treateed
        eigenvalues (armada vector&)  : Vector containing eigenvalues
        eigenvectors (armada matrix&) : Matrix containing eigenvalues
        maxiter (const int)           : selects the max allowed number of iterations
        counter (int&)                : counts the number of iterations
        converged (bool&)             : true if the algorithm has converged, false if not
        n (int)                       : size of side of matrix
    */
    double eps = 1e-6;
    double sum = 0;
    int k;
    int l; 
    while (!converged)
    {
        counter++;
        double max = largest_off_element(A, k, l, n);
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

void write_eig_to_file(arma::vec V, arma::mat R, int n, std::string filename)
{
    /*
    Writes the results from jacobi_eigensolver to a file
    Args: 
        V (armada vector) : vector containing eigenvalues
        R (armada vector) : matrix containing eigenvectors
        n (int)           : size of vector and matrix
        filename (string) : name of the file to write results to
    */
    std::fstream file;
    file.open(filename+".txt", std::ios::out);

    R = R.t();
    for (int i=0; i < n; i++)
    {
        file << V(i) << '\n';
        for (int j=0; j < n; j++)
        {

            file << R(i,j) << '\n';
        }
        file << '\n';
    }
    file.close();
}

void write_lin_to_file(arma::vec x, int n, std::string filename)
{
    /*
    Writes the linear vector x to a file
    Args:
        x (armada vector) : array to write to file
        n (int)           : size of vector
        filename (string) : name of file to write results to
    */
    std::fstream file;
    file.open(filename+".txt", std::ios::out);
    for (int i=0; i < n; i++)
    {
        file << x(i) << '\n';
    }
    file.close();
}
