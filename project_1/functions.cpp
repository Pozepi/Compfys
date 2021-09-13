#include "functions.hpp"



void linspace(double dummy[], double start, double end, int n)
{   
    /* 
    Fill dummy with linearly spaced values from start to end 
    with n data points 
    Args:
        dummy (double array): array to be filled with linspace values
        start (double)      : start value of the array
        end   (double)      : end value of the array
        n     (int)         : length of dummy 
    */
    dummy[0] = start;
    dummy[n] = end;
    double h = (end - start)/(n-1);
    for(int i = 1; i < n; i++)
    {
        dummy[i] = h + dummy[i-1];
    }
}

void Poisson_analytical(double dummy[], double x[], int n)
{   
    /*
    Create an array containing the analytical solution to the Poisson equation
    Args: 
        dummy (double array): array to be filled
        x     (double array): array containing x values used to calculate Poisson
        n     (int)         : length of both dummy and x
    */
    for (int i=0; i < n; i++)
    {
        dummy[i] = 1 - (1-expf(-10)) * x[i] - expf(-10 * x[i]);
    }
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

void fill_array(double dummy[], int n, double element)
{   
    /*
    Fills an input array with given value
    Args:
        dummy   (double array): array to be filled
        n       (int)         : length of dummy
        element (double )     : value of which to fill dummy with
    */
    for(int i = 0; i < n; i++)
    {
        dummy[i] = element;
    }
}

void solving_matrix(double dummy[], double x[], double a_array[], double b_array[], double c_array[], int n)
{
    /*
    Solvinges the Poisson eqiation numerically using the diagonals. Loops to fill out beta and gamma arrays which 
    are used to fill dummy with numerical Poisson values by looping backwards
    Args: 
        dummy   (double array) : array to be filled with numerical values
        x       (double array) : array containing x-values used to calculate numerical Poisson values
        a_array (double array) : subdiagonal array containing some values
        b_array (double array) : main diagonal 
        c_array (double array) : superdiagonal
        n       (int)          : length of dummy, x, a_array, b_array and c_array

    */
    double b = b_array[0];
    double beta[n];
    double gamma[n];
    double h = x[1]-x[0];

    beta[0] = b_array[0];
    gamma[0] = 0;

    // set up b~ and g
    for(int i = 1; i < n; i++)
    {
        beta[i]  = b_array[i] - (a_array[i]/beta[i-1])*c_array[i-1];
        gamma[i] = 100*exp(-10*x[i])*h*h - (a_array[i]/beta[i-1])*gamma[i-1];
    }    

    // fill out v_zeros
    dummy[n] = 0;
    dummy[0] = 0;
    
    for(int i = n-2; i > 0; i--)
    {
        dummy[i] = (gamma[i] - dummy[i+1]*c_array[i])/beta[i];
    }
}

void solving_special(double dummy[], double x[], int n)
{   
    /* 
    Special algorithm for solving the Poisson equation numerically with subdiagonal containing only -1, 
    main diagonal containing only 2 and superdiagonal containing only -1 with length n. Calculates the 
    numerical values by looping from 1 to n and fiding beta and gamma arrays then looping backwards
    to fill dummy.
    Args: 
        dummy (double array) : array to be filled with Poisson numerical values defined by x
        x     (double array) : array containing x values used to calculate Poisson numerical values
        n     (int)          : length of dummy and x
    */
    double a = -1;
    double b =  2;
    double c = -1;

    double h = x[1] - x[0];

    double beta[n];
    double gamma[n];

    beta[0]  = b;
    gamma[0] = 0;

    double ac_ = a*c;
    double hh_ = 100*h*h;

    for (int i = 1; i < n; i++)
    {
        beta[i]  = b - (ac_/beta[i-1]);
        gamma[i] = exp(-10*x[i])*hh_ - (a/beta[i-1])*gamma[i-1];
    }

    dummy[n] = 0;
    dummy[0] = 0;
    
    for(int j = n-1; j > 1; j--)
    {
        dummy[j] = (gamma[j] - dummy[j+1]*c)/beta[j];
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

void Delta(double dummy[], double u[], double v[], int n)
{   
    /*
    Calculates the absolute error between u and v of length n and fills it in dummy
    Args: 
        dummy (double array) : array to be filled with absolute error
        u     (double array) : first array used to calulate error 
        v     (double array) : second array used to calculate error
        n     (int)          : length of dummy, u and v
    */
    for (int i=0; i<n; i++)
    {
        dummy[i] = abs(u[i]-v[i]);
    }
}

void epsilon(double dummy[], double u[], double v[], int n)
{
    /* 
    Calculates the relative error between u and v, u being the true solution (in our case the
    analytical solution).
    Args: 
        dummy (double array) : array to be filled with relative error between u and v
        u     (double array) : array containing the "true" solution
        v     (doulbe array) : array used to calculate error by comparing to u
    */
    for (int i=0; i<n; i++)
    {
        dummy[i] = abs((u[i]-v[i])/u[i]);
    }
}
