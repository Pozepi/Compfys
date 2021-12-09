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
    int pot_type;
    double v0;
    std::string name;

    arma::cx_vec u; // internal matrix
    arma::cx_mat V;
    arma::cx_mat U;
    arma::sp_cx_mat A;
    arma::sp_cx_mat B;
    //int L_;
    //int N_;

    
    public:
    // Constructor
    Particle(double, double, double, double, double, double, double, double, double, int, double, std::string name_ = "");
    // Index transformer
    int transform_index(int, int); 
    // Make A and B vectors
    std::tuple<arma::sp_cx_mat, arma::sp_cx_mat> construct_AB();
    // update the system by one timestep
    void update_system();
    // perform simluation by running update_system() multiple times.
    void simulate_system();
    // set up initial state of u
    void initial_state();
    void potential(int slits);

    //arma::mat Create_lattice();
};


#endif