#include "penning_trap.cpp"

int main()
{
    PenningTrap penning_test;

    Particle one(1,2, {1,2,3}, {4,5,6});

    penning_test.add_particle(one);
    penning_test.add_n_particles(10, 1,2, {1,2,3}, {4,5,6});

    std::cout << "Number of particles " << penning_test.particle_count() << std::endl; 

    return 0;
}