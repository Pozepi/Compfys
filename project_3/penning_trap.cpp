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
// Force on particle i from particle j
arma::vec PenningTrap::force_particle(int i, int j)
{
    arma::vec F = {0,0,0};
    return F;
}
// Total force on particle i from external fields
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
// Total force on particle i from other particles
arma::vec PenningTrap::total_force_particles(int i)
{
    arma::vec F = {0,0,0};
    return F;
}
// Total force on particle i from both external fields and other particles
arma::vec PenningTrap::total_force(int i)
{
    double Fx = PenningTrap::total_force_external(i)(0)+PenningTrap::total_force_particles(i)(0);
    double Fy = PenningTrap::total_force_external(i)(1)+PenningTrap::total_force_particles(i)(1);
    double Fz = PenningTrap::total_force_external(i)(2)+PenningTrap::total_force_particles(i)(2);
    arma::vec F = {Fx,Fy,Fz};
    return F;
}

void PenningTrap::evolve_RK4(double dt)
{
    int null = 0;
}

void PenningTrap::evolve_forward_Euler(double dt, double time_stop)
{
    int size = PenningTrap::particle_count()*3;
    int N = time_stop/dt;
    arma::mat v(N, size);
    arma::mat pos(N, size);
    arma::vec time(N);
    
    int jump = 0;
    // setting initial values
    time(0) = 0;
    for(int k=0; k<PenningTrap::particle_count(); k++)
    {
        Particle particle = particles[k];
        v(0,jump) = particle.velocity()(0);
        v(0,jump+1) = particle.velocity()(1);
        v(0,jump+2) = particle.velocity()(2);
        pos(0,jump) = particle.position()(0);
        pos(0,jump+1) = particle.position()(1);
        pos(0,jump+2) = particle.position()(2);
        jump += 3;
    }

    for(int i=0; i<N-1; i++)
    {
        jump = 0;
        for(int j=0; j<PenningTrap::particle_count(); j++)
        {
            Particle particle_j = particles[j];

            arma::vec F = PenningTrap::total_force(j);
            arma::vec a = {F(0)/particle_j.mass(), F(1)/particle_j.mass(), F(2)/particle_j.mass()};
            v(i+1,jump) = v(i,jump) + a(0)*dt;
            v(i+1,jump+1) = v(i,jump+1) + a(1)*dt;
            v(i+1,jump+2) = v(i,jump+2) + a(2)*dt;
            pos(i+1,jump) = pos(i,jump)+v(i+1,jump)*dt;
            pos(i+1,jump+1) = pos(i,jump+1)+v(i+1,jump+1)*dt;
            pos(i+1,jump+2) = pos(i,jump+2)+v(i+1,jump+2)*dt;
            
            jump += 3;
        }
        // updating particles position and velocity for next iteration
        jump = 0;
        for(int k=0; k<PenningTrap::particle_count(); k++)
        {
            Particle particle = particles[0];
            double mass = particle.mass();
            double charge = particle.charge();

            arma::vec update_vel = {v(i+1,jump), v(i+1,jump+1), v(i+1,jump+2)};
            arma::vec update_pos = {pos(i+1,jump), pos(i+1,jump+1), pos(i+1,jump+2)};
            Particle update_particle(charge, mass, update_pos, update_vel);
            particles.erase(particles.begin());
            PenningTrap::add_particle(update_particle);

            jump += 3;
        }
        time(i+1) = time(i)+dt;
    }
}