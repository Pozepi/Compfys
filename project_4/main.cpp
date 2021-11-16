#include <iostream>
#include "lattice.hpp"

int main()
{
    // INITIALIZE PROGRAM
    std::cout << "Welcome to our program! Here are your options:" << '\n';
    std::cout << "  [1.] Simple lattice example (L=2, T=1)" << '\n';
    std::cout << "  [2.] Calculate Cv(T) and chi(T) for L=2" << '\n';
    
    std::cout << '\n' << "Please select an option" << "\n";
    int x;
    std::cin >> x;
    switch(x) 
    {
        case 1:
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

            break;
        }
        // loop over temperature
        // define my variables
        case 2: 
        {
            int L = 2;
            int N = L*L;
            double t1 = clock();
            double T0 = 0.001; double T1 = 100;
            int n = 1000; L = 2; N = L*L;
            double cycles = 200000;
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
                arma::vec average = myinstance.full_cycle(cycles);
                double Cv = myinstance.specific_heat_capacity(average);
                double chi = myinstance.susceptibility(average);

                vec_Cv(i) = Cv;
                vec_chi(i) = chi;

                std::cout << i << "/" << n << '\n';
            }

            //std::cout << vec_Cv << '\n';
            //std::cout << vec_chi << '\n';
            double t2 = clock();
            double duration_seconds = ((double) (t2-t1))/CLOCKS_PER_SEC;
            std::cout << "Time: " << duration_seconds << '\n';
            Temp.save("Temp.txt");
            vec_Cv.save("Cv.txt");
            vec_chi.save("chi.txt");

            break;
        }
        
        case 3:
            int noll = 0;
            break;
    }
    

    return 0;
}



