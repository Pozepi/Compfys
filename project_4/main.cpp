#include <iostream>
#include "lattice.hpp"

int main()
{
    int L = 2;
    int N = L*L;
    Lattice test(L, 1);
    test.Fill_lattice();
    int cycles = 10000;
    arma::vec average(6);
    average(5) = 0;

    clock_t t1 = clock();
    for(int i = 0; i < cycles; i++)
    {
        std::srand(time(0)+i);
        test.one_cycle_MCMC(N, average);
    }

    clock_t t2 = clock();
    double Cv = test.specific_heat_capacity(average);
    double chi = test.susceptibility(average);
    std::cout << "Cv: "<< Cv << '\n';
    std::cout << "chi: "<< chi << '\n';
    double duration_seconds = ((double) (t2-t1))/CLOCKS_PER_SEC;
    std::cout << "Time: " << duration_seconds << '\n';


    return 0;
}



