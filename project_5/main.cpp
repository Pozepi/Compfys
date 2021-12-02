#include <iostream>
#include "particle.hpp"

int main()
{
    int M = 100; // Remember inner size = (M-2, M-2)

    arma::vec a((M-2)*(M-2));
    a.ones();
    arma::vec b((M-2)*(M-2));
    b.ones();
    double h = 0.005;
    double dt = 0.1;
    // int M_, double h_, double dt_, double T_, 
    // double xc_, double yc_, double sigmax_, double sigmay_, double px_, double py_, 
    // arma::cx_vec v0_
    Particle myinstance(M, h, dt, 1., 0.5, 0.5, 0.5, 0.5, 0.1, 0.1, 1e5);
    myinstance.update_system();
    myinstance.potential(2);

    return 0;
}



