#include "penning_trap.hpp"

int main()
{
    double B0 = 9.65*10;
    double V0 = 9.65e+08;
    double d = 1e+04;
    double ke = 1.38935333e+05;

    PenningTrap penning_test(B0, V0, d);

    Particle one(1, 1, {1,2,3}, {4,5,6});
    Particle two(1,1, {0,2,1}, {2,1,7});
    Particle three(1,1,{5,8,9}, {5,1,2});

    penning_test.add_particle(one);
    penning_test.add_particle(two);
    penning_test.add_particle(three);

    penning_test.evolve_forward_Euler(0.1, 1, true);

    return 0;
}