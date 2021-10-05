#include "functions.hpp"
#include "particle.cpp"

class PenningTrap
{
    private:
    double magnetic_field_;
    double potential_;
    double dimension_;
    std::vector<Particle> particles;

    public:
    PenningTrap() {}

    PenningTrap(std::vector<Particle> particles_in)
    {
        int null = 0;
    }

    void add_particle(Particle particle_in)
    {
        particles.push_back(particle_in);
    }

    void add_n_particles(int n, double charge, double mass, arma::vec position, arma::vec velocity)
    {
        for(int i=0; i<n; i++)
        {
            Particle new_particle(charge, mass, position, velocity);
            add_particle(new_particle);
        }
    }

    double particle_count()
    {
        return particles.size();
    }

    void external_electric_field()
    {
        int null = 0;
    }

    void external_magnetic_field()
    {
        int null = 0;
    }

    void force_between_particles()
    {
        int null = 0;
    }

};