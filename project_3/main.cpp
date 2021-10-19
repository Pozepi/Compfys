#include "penning_trap.hpp"

int main()
{
    double volt = 9.64852558e+7;
    double cm = 1e+4;
    double B0 = 9.65*10;
    double V0 = 9.65e+08;
    double d = 1e+4;
    double ke = 1.38935333e+05;
    double charge = 1;
    double mass = 38.97;

    /*
    Particle one(charge, mass, {0,1,0}, {-1,0,0});
    Particle two(charge, mass, {0,-1,0}, {0,1,0});

    Particle one_z(charge, mass, {0,1,0}, {-1,0,1});
    Particle two_z(charge, mass, {0,-1,0}, {0,1,-1});

    PenningTrap penning_test_one1(B0, V0, d);
    penning_test_one1.add_particle(one);

    PenningTrap penning_test_one2(B0, V0, d);
    penning_test_one2.add_particle(one);

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

    PenningTrap penning_test_z1(B0, V0, d);
    penning_test_z1.add_particle(one_z);
    penning_test_z1.add_particle(two_z);

    PenningTrap penning_test_z2(B0, V0, d);
    penning_test_z2.add_particle(one_z);
    penning_test_z2.add_particle(two_z);

    penning_test_one1.evolve_forward_Euler(1e-3,100,false,true,"euler_one_particle");
    penning_test_one2.evolve_RK4(1e-3,100,false,true, "RK4_one_particle");
    penning_test.evolve_forward_Euler(1e-3, 100, true, true, "euler_two_particles_interaction");
    penning_test2.evolve_RK4(1e-3, 100, true, true, "RK4_two_particles_interaction");
    penning_test3.evolve_forward_Euler(1e-3, 100, false, true, "euler_two_particles_no_interaction");
    penning_test4.evolve_RK4(1e-3, 100, false, true, "RK4_two_particles_no_interaction");
    penning_test_z1.evolve_RK4(1e-3, 100, true, true, "RK4_two_particles_interaction_z_start");
    penning_test_z2.evolve_RK4(1e-3,100,false,true,"RK4_two_particles_no_interaction_z_start");
    */
    /*
    PenningTrap penning_test_against_analytical(B0, V0, d);

    Particle against_analytical(charge, mass, {1,0,1}, {0,1,0});
    penning_test_against_analytical.add_particle(against_analytical);
    penning_test_against_analytical.evolve_RK4(1e-3, 100, false, 0, 0, true, "RK4_against_analytical");
    */
    double d_10 = 0.05*cm;
    double V0_10 = 0.0025*volt;
    // double V0_time = V0*(1+f*std::cos(w_v*t));

    //PenningTrap penning_test_10(B0, V0_10, d_10);
    
    PenningTrap random_test(B0, V0_10, d_10);
    random_test.add_n_random_particles(10, charge, mass);
    random_test.evolve_RK4(1e-3, 100, true, 0,0,true, "RK4_random_test");
    std::cout << "Particles left in the penningtrap " << random_test.particles_inside_trap_count() << "\n";
    

    return 0;
}