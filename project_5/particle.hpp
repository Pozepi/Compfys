#ifndef __particle_hpp__  
#define __particle_hpp__
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
#include <complex>

class Particle
{
    private:
    // parameters of our instance
    int M; // matrix size with boundaru
    double h;
    double dt;
    double T; 
    double xc;
    double yc; 
    double sigmax; 
    double sigmay;
    double px;
    double py;
    arma::cx_mat v0; // idunno 
    arma::cx_vec u; // internal matrix
    arma::cx_mat V;
    arma::cx_mat A;
    arma::cx_mat B;
    //int L_;
    //int N_;

    
    public:
    // Constructor
    Particle(int, double, double, double, double, double, double, double, double, double, 
    arma::cx_vec);
    // Index transformer
    int transform_index(int, int); 
    // Make A and B vectors
    std::tuple<arma::cx_mat, arma::cx_mat> construct_AB();
    // update the system by one timestep
    void update_system();
    // set up initial state of u
    void initial_state();

    //arma::mat Create_lattice();
};


#endif