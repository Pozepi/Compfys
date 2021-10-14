#include "particle.hpp"
#include "penning_trap.hpp"


PenningTrap::PenningTrap(double magnetic_field, double potential, double dimension)
{
    magnetic_field_ = magnetic_field;
    potential_ = potential;
    dimension_ = dimension;
}

void PenningTrap::add_particle(Particle particle_in)
{
    particles.push_back(particle_in);
}

void PenningTrap::add_n_particles(int n, double charge, double mass, arma::vec position, arma::vec velocity)
{
    for(int i=0; i<n; i++)
    {
        Particle new_particle(charge, mass, position, velocity);
        add_particle(new_particle);
    }
}

double PenningTrap::particle_count()
{
    return particles.size();
}

arma::vec PenningTrap::external_E_field(arma::vec r)
{
    double x_der = 2*r[0]*potential_/(2*dimension_*dimension_);
    double y_der = 2*r[1]*potential_/(2*dimension_*dimension_);
    double z_der = -4*r[2]*potential_/(2*dimension_*dimension_);
    arma::vec E = {x_der,y_der,z_der};
    return E;
}

arma::vec PenningTrap::external_B_field(arma::vec r)
{
    arma::vec B = {0,0,magnetic_field_};
    return B;
}

arma::vec PenningTrap::force_particle(int i, int j)
{
    arma::vec F = {1,2,3};
    return F;
}

arma::vec PenningTrap::total_force_external(int i)
{
    Particle particle_i = particles[i];
    arma::vec v = particle_i.velocity();
    arma::vec B = PenningTrap::external_B_field(particle_i.position());
    arma::vec E = PenningTrap::external_E_field(particle_i.position());

    double Fx = particle_i.charge()*E(0)+particle_i.charge()*(v(1)*B(2)-v(2)*B(1));
    double Fy = particle_i.charge()*E(1)+particle_i.charge()*(v(0)*B(2)-v(2)*B(0));
    double Fz = particle_i.charge()*E(2)+particle_i.charge()*(v(0)*B(1)-v(1)*B(0));
    arma::vec F = {Fx, Fy, Fz};
    return F;
}

arma::vec PenningTrap::total_force_particles(int i)
{
    arma::vec F = {1,2,3};
    return F;
}

arma::vec PenningTrap::total_force(int i)
{
    arma::vec F = {1,2,3};
    return F;
}

void PenningTrap::evolve_RK4(double dt)
{
    int null = 0;
}

void PenningTrap::evolve_forward_Euler(double dt)
{
    int null = 0;
}