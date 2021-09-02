#include <iostream>
#include <math.h>
#include <fstream>
#include <iterator>

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

double * fill_array(double array[], int n, double element)
{
    for(int i = 0; i < n; i++)
    {
        array[i] = element;
    }
    return array;
}

double * Thomas(double lin[], double a_array[], double b_array[], double c_array[], double v_zeros[], int n)
{
    double b = b_array[0];
    double beta[n];
    double gamma[n];
    double h = lin[1]-lin[0];
    // std::cout << h;

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
        std::cout << i <<'\n';
    }
    
    return v_zeros;
}

void print_array(double array[], int n)
{
    for(int i = 0; i < n; i++)
    {
        std::cout << array[i] << '\n';
    }
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

int main()
{   
    int n;
    double a_value;
    double b_value;
    double c_value;

    std::cin >> n;
    std::cin >> a_value;
    std::cin >> b_value;
    std::cin >> c_value;

    double lin[n];
    double u_zeros[n];
    double a_zeros[n-1];
    double b_zeros[n];
    double c_zeros[n];
    double v_zeros[n];

    double * lins = linspace(0, 1, lin, n);
    std::cout << "-----linspace------" << '\n';
    print_array(lins, n);

    double * a_array = fill_array(a_zeros, n, a_value);
    double * b_array = fill_array(b_zeros, n, b_value);
    double * c_array = fill_array(c_zeros, n, c_value);

    std::cout << "-----Array print-----" << '\n';
    print_array(a_array, n);
    print_array(b_array, n);
    print_array(c_array, n);

    double * v_filled = Thomas(lins, a_array, b_array, c_array, v_zeros, n);

    
    std::cout << "-----V print-----" << '\n';
    print_array(v_filled, n);

    write_to_file(v_filled, lins, n);




    return 0;
}