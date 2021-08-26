#include <iostream>
#include <math.h>
#include <fstream>
#include <iterator>

double setting_array_elements (double a, double b, int N) // From a to b in N steps
{   
    double array[N];
    for (int i=0; i<N; i++)
        array[i] = a + (b-a)*i/N;
    return array;
}

double Poisson_analytical(double x[])
{   
    double u[sizeof(x)];
    for (int i=0; i < sizeof(x); i++)
    {
        u[i] = 1 - (1-exp(-10)) * x[i] - exp(-10 * x[i]);
    }
    return u;
}

void write_to_file(double a1[], double a2[])
{
    std::fstream file;
    file.open("values.txt", std::ios::out);
    for (int i=0; i<sizeof(a1); i++)
    {
        file << a1[i] << ' ' << a2[i] << '\n';
    }
    file.close();
}

int main()
{
    double x = setting_array_elements(0, 30, 10);
    double u = Poisson_analytical(x);
    write_to_file(x, u);
    return 0;
}
