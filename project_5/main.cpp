#include <iostream>
#include "particle.hpp"

int main()
{
    double h = 0.005;
    double dt = 2.5e-5;
    // int M_, double h_, double dt_, double T_, 
    // double xc_, double yc_, double sigmax_, double sigmay_, double px_, double py_, 
    // arma::cx_vec v0_

    // h, dt, T, xc, yc, sigmax, sigmay, px, py, pot_type, v0, name

    // double slit with no potential
    // Particle myinstance(h, dt, 0.008, 0.25, 0.5, 0.05, 0.05, 200, 0.0, 2, 0, "7");
    Particle myinstance(h, dt, 0.002, 0.25, 0.5, 0.05, 0.1, 200, 0.0, 2, 1e10, "7_2");

    // double slit with potential
    //Particle myinstance(h, dt, 0.002, 0.25, 0.5, 0.05, 0.2, 200, 0.0, 2, 1e10, "8");

    // no slit
    // Particle myinstance(h, dt, 0.008, 0.25, 0.5, 0.05, 0.05, 200, 0.0, 0, 1e4);
    myinstance.simulate_system();

    return 0;
}



