#include <iostream>
#include <math.h>
#include <armadillo>


float Poisson(arma::rowvec x, arma::rowvec u, int size)
{
    for (int i=0; i<size; i++)
    {
        u[i] = 1 - (1-exp(-10)) * x[i] - exp(-10 * x[i]);
        std::cout << u[i];
    }
    return u
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

    arma::rowvec x = arma::linspace<arma::rowvec>(start,end,size);
    arma::rowvec u = arma::zeros<arma::rowvec>(size);
    Poisson(x, size);
    std::cout << x;

}