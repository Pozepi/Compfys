#ifndef __lattice_hpp__  
#define __lattice_hpp__
#include <armadillo>
#include <iostream>
#include <sys/stat.h>
#include <string>
#include <cstdlib>
#include <time.h>
#include <numeric>
#include <random>
#include <algorithm>
#include <map>
#include <time.h>
#include "omp.h"
#include <chrono>

class Lattice
{
    private:
    int L_;
    int N_;
    double T_;
    double E;
    double M;
    double dE;
    double r;
    double p;

    arma::imat lattice;
    arma::mat padded;
    std::map<double, double> my_map;

    
    public:
    Lattice(int L, double T, bool ordered);
    arma::mat Create_lattice();
    void Fill_lattice(bool ordered);

    arma::mat Pad_lattice(arma::mat lat);
    arma::mat unpad(arma::mat pad);

    double Total_magnetization();
    double Total_energy();
    double Boltzman();

    arma::mat Replace_pad(arma::mat lat);

    double energy_per_spin();
    double energy_per_spin_expectation(arma::vec average);

    double magnetization_per_spin();
    double magnetization_per_spin_expectation(arma::vec average);

    double specific_heat_capacity(arma::vec average);
    double susceptibility(arma::vec average);
    void one_cycle_MCMC(arma::vec& average);
    void burn_in(int cycles);
    bool test_flip(int i,int j);
    int periodic(int i, int limit, int add);
    // void one_cycle_MCMC_2(arma::vec& average, std::map<double, double> my_map);

    // double calc_dE_2(arma::mat S1, arma::mat S2, int i, int j);
    // double calc_dE_N(arma::mat S1, arma::mat S2, int i, int j);
    arma::vec full_cycle(int cycles, arma::vec& eps_list, arma::vec& m_list, bool sample_eps_lattice=false);
};


#endif