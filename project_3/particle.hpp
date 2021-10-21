#ifndef __particle_hpp__  
#define __particle_hpp__

#include <armadillo>
#include <iostream>
#include <sys/stat.h>
#include <string>
#include <cstdlib>


class Particle
{
    private:

    double charge_;
    double mass_;
    arma::vec position_;
    arma::vec velocity_;

    public:

    Particle(double charge, double mass, arma::vec position, arma::vec velocity);
    double charge();
    double mass();
    arma::vec position();
    arma::vec velocity();
};

#endif