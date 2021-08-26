#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <armadillo>

float * Poisson_analytical(arma::vec x_array, float u_array[], int size)
{
    for (int i=0; i<size; i++)
    {
        u_array[i] = 1 - (1-exp(-10)) * x_array[i] - exp(-10 * x_array[i]);
        std::cout << std::setprecision(6) << x_array[i] << ' ' << u_array[i];
        std::cout << '\n';
    }
    return u_array;
}

void write_to_file(arma::vec array, float poisson[], int size)
{
    std::fstream file;
    file.open("values.txt", std::ios::out);
    for (int i=0; i<size; i++)
    {
        file << std::setprecision(6) << array[i] << ' ' << poisson[i] << '\n';
    }
    file.close();
}

int main()
{
    int start;
    int end;
    int size;

    std::cout << "Start value: ";
    std::cin >> start;
    std::cout << "End value: ";
    std::cin >> end;
    std::cout << "Array size: ";
    std::cin >> size;

    arma::vec x_array = arma::linspace(start, end, size);
    float u_array[size];

    float * u = Poisson_analytical(x_array, u_array, size);
    write_to_file(x_array, u, size);

    return 0;
}
