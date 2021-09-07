#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>

float * linspace(float a, float b, int N) // From a to b in N steps
{   
    float array[N];
    for (int i=0; i<N; i++)
        array[i] = a + (b-a)*i/N;
    return array;
}


float * Poisson_analytical(float x[])
{   
    float u[sizeof(x)];
    for (int i=0; i<sizeof(x); i++)
    {
        u[i] = 1 - (1-exp(-10)) * x[i] - exp(-10 * x[i]);
        //std::cout << std::setprecision(6) << x[i] << ' ' << u[i];
        //std::cout << '\n';
    }
    return u;
}

void write_to_file(float array[], float poisson[])
{
    std::fstream file;
    file.open("values.txt", std::ios::out);
    for (int i=0; i<sizeof(array); i++)
    {
        file << std::setprecision(6) << array[i] << ' ' << poisson[i] << '\n';
    }
    file.close();
}

int main()
{
    float start;
    float end;
    int length;

    float x[length];
    float u[length];

    x = linspace(start, end, length);

    u = Poisson_analytical(x);
    write_to_file(x, u);

    

    return 0;
}
