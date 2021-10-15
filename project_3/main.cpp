#include "penning_trap.hpp"

int main()
{
    double B0 = 9.65*10;
    double V0 = 9.65e+08;
    double d = 1e+04;
    double ke = 1.38935333e+05;
    double charge = 1;
    double mass = 38.97;

    Particle one(charge, mass, {0,1,0}, {-1,0,0});
    Particle two(charge, mass, {0,-1,0}, {0,1,0});

    Particle one_z(charge, mass, {0,1,0}, {-1,0,1});
    Particle two_z(charge, mass, {0,-1,0}, {0,1,-1});

    PenningTrap penning_test_one1(B0, V0, d);
    penning_test_one1.add_particle(one);

    PenningTrap penning_test_one2(B0, V0, d);
    penning_test_one2.add_particle(one)

    PenningTrap penning_test(B0, V0, d);
    penning_test.add_particle(one);
    penning_test.add_particle(two);

    PenningTrap penning_test2(B0, V0, d);
    penning_test2.add_particle(one);
    penning_test2.add_particle(two);

    PenningTrap penning_test3(B0, V0, d);
    penning_test3.add_particle(one);
    penning_test3.add_particle(two);

    PenningTrap penning_test4(B0, V0, d);
    penning_test4.add_particle(one);
    penning_test4.add_particle(two);

    penning_test_one1.evolve_forward_Euler(1e-3,100,false,true,"euler_one_particle");
    penning_test_one2.evolve_RK4(1e-3,100,false,true )
    penning_test.evolve_forward_Euler(1e-3, 100, true, true, "euler_two_particles_interaction");
    penning_test2.evolve_RK4(1e-3, 100, true, true, "RK4_two_particles_interaction");
    penning_test3.evolve_forward_Euler(1e-3, 100, false, true, "euler_two_particles_no_interaction");
    penning_test4.evolve_RK4(1e-3, 100, false, true, "RK4_two_particles_interaction");


    return 0;
}