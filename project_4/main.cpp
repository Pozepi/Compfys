#include <iostream>
#include "lattice.hpp"

int main()
{
    int L = 2;
    int N = L*L;
    Lattice test(L, 1);
    test.Fill_lattice();

    clock_t t1 = clock();
    arma::vec average = test.full_cycle(100000);
    clock_t t2 = clock();
    double Cv = test.specific_heat_capacity(average);
    double chi = test.susceptibility(average);
    std::cout << "Cv: "<< Cv << '\n';
    std::cout << "chi: "<< chi << '\n';
    double duration_seconds = ((double) (t2-t1))/CLOCKS_PER_SEC;
    std::cout << "Time: " << duration_seconds << '\n';

    // loop over temperature
    // define my variables
    
    t1 = clock();
    double T0 = 0.001; double T1 = 100;
    int n = 100; L = 2; N = L*L;
    double cycles = 100000;
    // define my vectors
    //arma::vec Temp = arma::linspace(T0, T1, n);
    arma::vec Temp = arma::logspace(-2, 2, n);
    arma::vec vec_Cv(n);
    arma::vec vec_chi(n);
    // loop and fill vectors
    for (int i=0; i<n; i++)
    {
        double Tempi = Temp(i);
        Lattice myinstance(L, Tempi);
        average = myinstance.full_cycle(cycles);
        Cv = myinstance.specific_heat_capacity(average);
        chi = myinstance.susceptibility(average);

        vec_Cv(i) = Cv;
        vec_chi(i) = chi;

        std::cout << i << "/" << n << '\n';
    }

    //std::cout << vec_Cv << '\n';
    //std::cout << vec_chi << '\n';
    t2 = clock();
    duration_seconds = ((double) (t2-t1))/CLOCKS_PER_SEC;
    std::cout << "Time: " << duration_seconds << '\n';
    Temp.save("Temp.txt");
    vec_Cv.save("Cv.txt");
    vec_chi.save("chi.txt");
    

    return 0;
}



