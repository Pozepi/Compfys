#include <iostream>
#include <math.h>

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

double * Thomas(double u_array[], double a_array[], double b_array[], double c_array[], double x_zeros[], int n)
{
    double b = b_array[0];
    double alpha[n];
    double beta[n];

    alpha[0] = b_array[0];

    for(int i = 1; i < n; i++)
    {
        alpha[i] = b_array[i] - a_array[i-1]/alpha[i-1]*c_array[i-1];
    }    
    return 0;

}

void print_array(double array[], int n)
{
    for(int i = 0; i < n; i++)
    {
        std::cout << array[i] << '\n';
    }
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
    double x_zeros[n];

    double * lins = linspace(0, 1, lin, n);
    std::cout << "-----linspace------" << '\n';
    print_array(lins, n);

    double * u_array = Poisson_analytical(lins, u_zeros, n);
    double * a_array = fill_array(a_zeros, n-1, a_value);
    double * b_array = fill_array(b_zeros, n, b_value);
    double * c_array = fill_array(c_zeros, n, c_value);

    std::cout << "-----Array print-----" << '\n';
    print_array(a_array, n-1);
    print_array(b_array, n);
    print_array(c_array, n);



    return 0;
}