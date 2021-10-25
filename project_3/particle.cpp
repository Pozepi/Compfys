#include "particle.hpp"

Particle::Particle(double charge, double mass, arma::vec position, arma::vec velocity)
{
    /*
    CONSTRUCTOR

    Creates a particle with charge, mass, position and velocity.
    Args:
        charge      (double)    :   Charge of particle
        mass        (double)    :   mass of particle
        position    (arma::vec) :   initial position of particle (x,y,z)
        velocity    (arma::vec) :   initial velocity of particle (x,y,z)
    */
    charge_ = charge;
    mass_ = mass;
    position_ = position;
    velocity_ = velocity;
}

double Particle::charge()
{
    /*
    Particle charge
    Returns:
        charge_     (double)    :   Charge of particle
    */
    return charge_;
}

double Particle::mass()
{
    /*
    Particle mass
    Returns:
        mass_   (double)    :   mass of particle
    */
    return mass_;
}

arma::vec Particle::position()
{
    /*
    Particle position
    Returns:
        position_   (arma::vec) :   particle position in (x,y,z)
    */
    return position_;
}

arma::vec Particle::velocity()
{
    /*
    Particle velocity
    Returns:
        velocity_   (arma::vec) :   particle velocity in (x,y,z)
    */
    return velocity_;
}
