#include "particle.hpp"

Particle::Particle(double charge, double mass, arma::vec position, arma::vec velocity)
{
    charge_ = charge;
    mass_ = mass;
    position_ = position;
    velocity_ = velocity;
}

double Particle::charge()
{
    return charge_;
}

double Particle::mass()
{
    return mass_;
}

arma::vec Particle::position()
{
    return position_;
}

arma::vec Particle::velocity()
{
    return velocity_;
}
