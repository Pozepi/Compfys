#include "penning_trap.hpp"

int main()
{
    double B0 = 9.65*10;
    double V0 = 9.65e+08;
    double d = 1e+04;
    double ke = 1.38935333e+05;
    double charge = 1;
    double mass = 38.97;

    PenningTrap penning_test(B0, V0, d);
    PenningTrap penning_test2(B0, V0, d);

    Particle one(charge, mass, {0,1,0}, {-1,0,0});
    Particle two(1,1, {-1,0,0}, {0,1,0});
    //Particle three(1,1,{5,8,9}, {5,1,2});

    penning_test.add_particle(one);
    penning_test2.add_particle(one);
    //penning_test.add_particle(two);
    //penning_test2.add_particle(two);

    penning_test.evolve_forward_Euler(1e-3, 100, false, true, "euler_one_particle");

    penning_test2.evolve_RK4(1e-3, 100, false, true, "RK4_one_particle");

    return 0;
}