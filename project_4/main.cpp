#include <iostream>
#include "lattice.hpp"

int main()
{
    int L = 2;
    int N = L*L;
    Lattice test(L, 1);
    test.Fill_lattice();
    int cycles = 10000;
    double eps;
    double m;
    double Cv;
    double chi;

    arma::vec eps_list(cycles);
    arma::vec m_list(cycles);

    clock_t t1 = clock();
    for(int i = 0; i < cycles; i++)
    {
        std::srand(time(0)+i);
        test.one_cycle_MCMC(14, eps, m);
        eps_list(i) = eps;
        m_list(i) = m;
    }
    Cv = test.specific_heat_capacity(eps_list);
    chi = test.susceptibility(m_list);
    clock_t t2 = clock();
    std::cout << "CV: " << Cv << '\n';
    std::cout << "Chi: " << chi << '\n';
    eps_list.save("eps.txt");
    m_list.save("m.txt");
    double duration_seconds = ((double) (t2-t1))/CLOCKS_PER_SEC;
    std::cout << "Time: " << duration_seconds << '\n';


    return 0;
}



