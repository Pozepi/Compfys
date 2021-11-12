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
    double Total_magnetization();
    double Total_energy(arma::mat lat);
    double Boltzman(arma::mat lat);
    arma::mat Replace_pad(arma::mat lat);
    void one_cycle_MCMC(int n);
};


#endif