#include "functions.hpp"

// Each function requires to input a prebuilt
// array to be modified with output values.

void linspace(double dummy[], double start, double end, int n)
{   
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
    for (int i=0; i < n; i++)
    {
        dummy[i] = 1 - (1-expf(-10)) * x[i] - expf(-10 * x[i]);
    }
}

void print_array(double array[], int n)
{
    for(int i = 0; i < n; i++)
    {
        std::cout << array[i] << '\n';
    }
}

void fill_array(double dummy[], int n, double element)
{
    for(int i = 0; i < n; i++)
    {
        dummy[i] = element;
    }
}

void solving_matrix(double dummy[], double x[], double a_array[], double b_array[], double c_array[], int n)
{
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

void write_to_file(std::string filename, double array[], double lin[], int n)
{
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
    for (int i=0; i<n; i++)
    {
        dummy[i] = abs(u[i]-v[i]);
        //std::cout << std::setprecision(6) << x[i] << ' ' << u[i];
        //std::cout << '\n';
    }
}

void epsilon(double dummy[], double u[], double v[], int n)
{
    for (int i=1; i<n-1; i++)
    {
        dummy[i] = abs((u[i]-v[i])/u[i]);
    }
}
