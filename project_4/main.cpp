#include <iostream>
#include "lattice.hpp"

int main()
{
    int L = 2;
    int N = L*L;
    Lattice test(L, 1);
    test.Fill_lattice();

    clock_t t1 = clock();
    arma::vec average = test.full_cycle(1000);
    clock_t t2 = clock();
    double Cv = test.specific_heat_capacity(average);
    double chi = test.susceptibility(average);
    std::cout << "Cv: "<< Cv << '\n';
    std::cout << "chi: "<< chi << '\n';
    double duration_seconds = ((double) (t2-t1))/CLOCKS_PER_SEC;
    std::cout << "Time: " << duration_seconds << '\n';

    // loop over temperature
    // define my variables
    
    auto t3 = std::chrono::steady_clock::now();
    double T0 = 0.001; double T1 = 100;
    int n = 100; L = 10; N = L*L;
    double cycles = 10000;
    // define my vectors
    //arma::vec Temp = arma::linspace(T0, T1, n);
    arma::vec Temp = arma::logspace(-2, 2, n);
    arma::vec vec_Cv(n);
    arma::vec vec_chi(n);
    // loop and fill vectors
    #pragma omp parallel for 
    for (int i=0; i<n; i++)
    {
        double Tempi = Temp(i);
        Lattice myinstance(L, Tempi);
        average = myinstance.full_cycle(cycles);
        Cv = myinstance.specific_heat_capacity(average);
        chi = myinstance.susceptibility(average);
        vec_Cv(i) = Cv;
        vec_chi(i) = chi;
    //std::cout << i << "/" << n << '\n';
    }

    auto t4 = std::chrono::steady_clock::now();
    std::cout<<std::chrono::duration_cast<std::chrono::milliseconds>(t4-t3).count()<<" milliseconds"<<"\n";
    Temp.save("Temp.txt");
    vec_Cv.save("Cv.txt");
    vec_chi.save("chi.txt");
    

    return 0;
}



