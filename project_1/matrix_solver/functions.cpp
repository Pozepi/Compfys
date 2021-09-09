#include "functions.hpp"
#include <math.h>  // ??

double * linspace(double start, double end, double lin[], int n)
{
    lin[0] = start;
    lin[n] = end;
    double h = (end - start)/(n-1);
    for(int i = 1; i < n; i++)
    {
        lin[i] = h + lin[i-1];
    }
    return lin;
}

double * Poisson_analytical(double x[])
{   
    double u[sizeof(x)];
    for (int i=0; i<sizeof(x); i++)
    {
        u[i] = 1 - (1-exp(-10)) * x[i] - exp(-10 * x[i]);
        //std::cout << std::setprecision(6) << x[i] << ' ' << u[i];
        //std::cout << '\n';
    }
    return u;
}

void print_array(double array[], int n)
{
    for(int i = 0; i < n; i++)
    {
        std::cout << array[i] << '\n';
    }
}

double * fill_array(double array[], int n, double element)
{
    for(int i = 0; i < n; i++)
    {
        array[i] = element;
    }
    return array;
}

double * solving_matrix(double lin[], double a_array[], double b_array[], double c_array[], double v_zeros[], int n)
{
    double b = b_array[0];
    double beta[n];
    double gamma[n];
    double h = lin[1]-lin[0];

    beta[0] = b_array[0];
    gamma[0] = 0;

    // set up b~ and g
    for(int i = 1; i < n; i++)
    {
        beta[i]  = b_array[i] - (a_array[i]/beta[i-1])*c_array[i-1];
        gamma[i] = 100*exp(-10*lin[i])*h*h - (a_array[i]/beta[i-1])*gamma[i-1];
    }    

    // fill out v_zeros
    v_zeros[n] = 0;
    v_zeros[0] = 0;
    
    for(int i = n-2; i > 0; i--)
    {
        v_zeros[i] = (gamma[i] - v_zeros[i+1]*c_array[i])/beta[i];
    }
    
    return v_zeros;
}

void write_to_file(double array[], double lin[], int n)
{
    std::fstream file;
    file.open("values2.txt", std::ios::out);
    for (int i=0; i<n; i++)
    {
        file << lin[i] << ' ' << array[i] << '\n';
    }
    file.close();
}

double * Delta(double u[], double v[])
{   
    double error[sizeof(u)];
    for (int i=0; i<sizeof(u); i++)
    {
        error[i] = log(abs(u[i]-v[i]));
        //std::cout << std::setprecision(6) << x[i] << ' ' << u[i];
        //std::cout << '\n';
    }
    return error;
}

double * epsilon(double u[], double v[])
{
    double eps[sizeof(u)];
    for (int i=0; i<sizeof(u); i++)
    {
        eps[i] = log(abs((u[i]-v[i])/u[i]));
    }
    return eps;
}