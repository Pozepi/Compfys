#include <iostream>
#include "lattice.hpp"

int main()
{
    Lattice test(2, 1);
    test.Fill_lattice();
    clock_t t1 = clock();
    int cycles = 10000;
    double eps;
    double m;
    double Cv;
    double chi;

    arma::vec eps_list(cycles);
    arma::vec m_list(cycles);
    arma::vec Cv_list(cycles);
    arma::vec chi_list(cycles);


    for(int i = 0; i < cycles; i++)
    {
        std::srand(time(0)+i);
        test.one_cycle_MCMC(4, eps, m, Cv, chi);
        eps_list(i) = eps;
        m_list(i) = m;
        Cv_list(i) = Cv;
        chi_list(i) = chi;
    }
    clock_t t2 = clock();
    eps_list.save("eps.txt");
    m_list.save("m.txt");
    double duration_seconds = ((double) (t2-t1))/CLOCKS_PER_SEC;
    std::cout << duration_seconds << '\n';


    return 0;
}



