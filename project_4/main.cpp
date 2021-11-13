#include <iostream>
#include "lattice.hpp"

int main()
{
    int L = 2;
    int N = L*L;
    Lattice test(L, 1);
    test.Fill_lattice();
    clock_t t1 = clock();
    int cycles = 1000;
    double eps;
    double m;
    double Cv;
    double chi;

    arma::vec eps_list(cycles);
    arma::vec m_list(cycles);

    for(int i = 0; i < cycles; i++)
    {
        std::srand(time(0)+i);
        test.one_cycle_MCMC(std::pow(2, N), eps, m);
        eps_list(i) = eps;
        m_list(i) = m;
    }
    Cv = test.specific_heat_capacity(eps_list);
    std::cout << "CV: " << Cv << '\n';
    clock_t t2 = clock();
    eps_list.save("eps.txt");
    m_list.save("m.txt");
    double duration_seconds = ((double) (t2-t1))/CLOCKS_PER_SEC;
    std::cout << "Time: " << duration_seconds << '\n';


    return 0;
}



