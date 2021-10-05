#include "functions.hpp"

class Particle
{
    private:
    double charge_;
    double mass_;
    arma::vec position_;
    arma::vec velocity_;

    public:
    Particle(double charge, double mass, arma::vec position, arma::vec velocity)
    {
        charge_ = charge;
        mass_ = mass;
        position_ = position;
        velocity_ = velocity;
    }

    double charge()
    {
        return charge_;
    }

    double mass()
    {
        return mass_;
    }

    arma::vec position()
    {
        return position_;
    }

    arma::vec velocity()
    {
        return velocity_;
    }

};