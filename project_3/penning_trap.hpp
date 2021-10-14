#ifndef __penning_trap_hpp__  
#define __penning_trap_hpp__

#include "particle.hpp"

class PenningTrap
{
    private:

    double magnetic_field_;
    double potential_;
    double dimension_;
    std::vector<Particle> particles;

    public:

    PenningTrap() {}

    PenningTrap(double magnetic_field, double potential, double dimension);

    void add_particle(Particle particle_in);

    void add_n_particles(int n, double charge, double mass, arma::vec position, arma::vec velocity);

    double particle_count();
    
    arma::vec external_E_field(arma::vec r);

    arma::vec external_B_field(arma::vec r);

    arma::vec force_particle(int i, int j);

    arma::vec total_force_external(int i);

    arma::vec total_force_particles(int i); 

    arma::vec total_force(int i);

    void evolve_RK4(double dt, double time_stop, bool makefile=false, std::string filename="RK4");

    void evolve_forward_Euler(double dt, double time_stop, bool makefile=false, std::string filename="forward_euler");
};

#endif


