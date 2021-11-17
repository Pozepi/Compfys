#include <iostream>
#include "lattice.hpp"

void loop_over_temp(int L, std::string filename, int cycles)
{
    int N = L*L;
    // loop over temperature
    // define my variables
    
    auto t1 = std::chrono::steady_clock::now();
    double T0 = 2.1; double T1 = 2.4;
    int n = 100;
    // define my vectors
    arma::vec Temp = arma::logspace(-2, 2, n);
    arma::vec vec_Cv(n);
    arma::vec vec_chi(n);
    arma::vec vec_eps(n);
    arma::vec vec_m(n);
    // loop and fill vectors
    #pragma omp parallel for 
    for (int i=0; i<n; i++)
    {
        double Tempi = Temp(i);
        Lattice myinstance(L, Tempi);
        arma::vec average = myinstance.full_cycle(cycles);
        double Cv = myinstance.specific_heat_capacity(average);
        double chi = myinstance.susceptibility(average);
        double eps = myinstance.energy_per_spin_expectation(average);
        double m = myinstance.magnetization_per_spin_expectation(average);
        vec_Cv(i) = Cv;
        vec_chi(i) = chi;
        vec_eps(i) = eps;
        vec_m(i) = m;
    }

    auto t2 = std::chrono::steady_clock::now();
    std::cout<<std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()<<" milliseconds"<<"\n";
    Temp.save("Temp"+filename);
    vec_Cv.save("Cv"+filename);
    vec_chi.save("chi"+filename);
    vec_eps.save("expect_eps"+filename);
    vec_m.save("expect_m"+filename);
}

int main()
{
    // INITIALIZE PROGRAM
    std::cout << "Welcome to our program! Here are your options:" << '\n';
    std::cout << "  [1.] Simple lattice example (L=2, T=1)" << '\n';
    std::cout << "  [2.] Calculate Cv(T) and chi(T) for L=2" << '\n';
    std::cout << "  [3.] Calculate Cv(T) and chi(T) for L=20" << '\n';
    std::cout << "  [4.] Calculate Cv(T) and chi(T) for L=40" << '\n';
    std::cout << "  [5.] Calculate Cv(T) and chi(T) for L=60" << '\n';
    std::cout << "  [6.] Calculate Cv(T) and chi(T) for L=80" << '\n';
    std::cout << "  [7.] Calculate Cv(T) and chi(T) for L=100" << '\n';
    std::cout << "  [8.] Calculate Cv(T) and chi(T) for L=n and cycles=m" << '\n';

    
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
            double cycles = 50000;
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
        {
            loop_over_temp(20, "_L20.txt", 100000);
            break;
        }

        case 4: 
        {
            loop_over_temp(40, "_L40.txt", 100000);
            break;
        }

        case 5:
        {
            loop_over_temp(60, "_L60.txt", 100000);
            break;
        }

        case 6:
        {
            loop_over_temp(80, "_L80.txt", 100000);
            break;
        }

        case 7:
        {
            loop_over_temp(100, "_L100.txt", 100000);
            break;
        }

        case 8:
        {   
            int L;
            int cycles;
            std::cout << "Input lattice size: ";
            std::cin >> L;
            std::cout << "Input number of Monte Carlo cycles: ";
            std::cin >> cycles;
            loop_over_temp(L, "_Ln.txt", cycles);
            break;
        }
        
    }


    return 0;
}



