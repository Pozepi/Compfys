#include "penning_trap.hpp"

int main()
{
    double B0 = 9.65*10;
    double V0 = 9.65e+08;
    double d = 1e+04;
    double ke = 1.38935333e+05;

    PenningTrap penning_test(B0, V0, d);

    Particle one(1, 2, {1,2,3}, {4,5,6});

    penning_test.add_particle(one);

    std::cout << "Number of particles " << penning_test.particle_count() << std::endl;
    std::cout << "Force on particle: \n" << penning_test.total_force_external(0);

    return 0;
}