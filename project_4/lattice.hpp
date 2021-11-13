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

class Lattice
{
    private:
    int L_;
    int N_;
    double T_;
    arma::mat lattice;
    arma::mat padded;
    
    public:
    Lattice(int L, double T);
    arma::mat Create_lattice();
    void Fill_lattice();
    arma::mat Pad_lattice(arma::mat lat);
    arma::mat unpad(arma::mat pad);
    double Total_magnetization(arma::mat lat, bool padded);
    double Total_energy(arma::mat lat, bool padded);
    double Boltzman();
    arma::mat Replace_pad(arma::mat lat);
    double energy_per_spin(arma::mat lat, bool padded);
    double magnetization_per_spin(arma::mat lat, bool padded);
    double specific_heat_capacity(arma::vec eps);
    double susceptibility();
    void one_cycle_MCMC(int n, double& eps, double& m);
};


#endif