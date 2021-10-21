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


    arma::vec h = {1e-1, 1e-2, 1e-3, 1e-4, 1e-5};
    Particle one(charge, mass, {1,0,1}, {0,1,0});

    PenningTrap h_I(B0,V0,d);
    PenningTrap h_II(B0,V0,d);
    PenningTrap h_III(B0,V0,d);
    PenningTrap h_IV(B0,V0,d);
    PenningTrap h_V(B0,V0,d);
    PenningTrap h_VI(B0,V0,d);
    PenningTrap h_VII(B0,V0,d);
    PenningTrap h_VIII(B0,V0,d);
    PenningTrap h_IX(B0,V0,d);
    PenningTrap h_X(B0,V0,d);

    h_I.add_particle(one);
    h_II.add_particle(one);
    h_III.add_particle(one);
    h_IV.add_particle(one);
    h_V.add_particle(one);
    h_VI.add_particle(one);
    h_VII.add_particle(one);
    h_VIII.add_particle(one);
    h_IX.add_particle(one);
    h_X.add_particle(one);

    h_I.evolve_RK4(h(0),100,false,0,0,true,"an_RK4_he1");
    h_II.evolve_RK4(h(1),100,false,0,0,true,"an_RK4_he2");
    h_III.evolve_RK4(h(2),100,false,0,0,true,"an_RK4_he3");
    h_IV.evolve_RK4(h(3),100,false,0,0,true,"an_RK4_he4");
    h_V.evolve_RK4(h(4),100,false,0,0,true,"an_RK4_he5");

    h_VI.evolve_forward_Euler(h(0),100,false,0,0,true,"an_EU_he1");
    h_VII.evolve_forward_Euler(h(1),100,false,0,0,true,"an_EU_he2");
    h_VIII.evolve_forward_Euler(h(2),100,false,0,0,true,"an_EU_he3");
    h_IX.evolve_forward_Euler(h(3),100,false,0,0,true,"an_EU_he4");
    h_X.evolve_forward_Euler(h(4),100,false,0,0,true,"an_EU_he5");

    
    /*
    Particle nineone(charge, mass, {1,0,1}, {0,1,0});

    PenningTrap nouno(B0, V0, d);
    PenningTrap nodos(B0, V0, d);
    PenningTrap notres(B0, V0, d);
    PenningTrap noquatro(B0, V0, d);

    nouno.add_particle(nineone);
    nodos.add_particle(nineone);
    notres.add_particle(nineone);
    noquatro.add_particle(nineone);

    nouno.evolve_RK4(1e-3, 100, false, 0, 0, true, "nineone_RK4_te3");
    nodos.evolve_forward_Euler(1e-3, 100, false, 0, 0, true, "nineone_EU_te3");
    notres.evolve_RK4(1e-4, 100, false, 0, 0, true, "nineone_RK4_te4");
    noquatro.evolve_forward_Euler(1e-4,100,false, 0,0,true, "nineone_EU_te4");
    */
    /*
    Particle ninetwo(charge,mass,{1,0,0},{0,1,0});
    Particle ninetwotwo(charge,mass,{-1,0,0},{1,0,0});

    PenningTrap ntuno(B0, V0, d);
    PenningTrap ntdos(B0, V0, d);
    PenningTrap nttres(B0, V0, d);
    PenningTrap ntquatro(B0, V0, d);
    PenningTrap ntcinco(B0, V0, d);

    ntuno.add_particle(ninetwo);
    ntuno.add_particle(ninetwotwo);
    ntdos.add_particle(ninetwo);
    ntdos.add_particle(ninetwotwo);
    nttres.add_particle(ninetwo);
    nttres.add_particle(ninetwotwo);
    ntquatro.add_particle(ninetwo);
    ntquatro.add_particle(ninetwotwo);
    ntcinco.add_particle(ninetwo);
    ntcinco.add_particle(ninetwotwo);

    ntuno.evolve_RK4(1e-3, 100, false, 0, 0, true, "ninetwo_RK4_te3");
    ntdos.evolve_RK4(1e-4, 100, false, 0, 0, true, "ninetwo_RK4_te4");
    nttres.evolve_forward_Euler(1e-3, 100, false, 0, 0, true, "ninetwo_EU_te3");
    ntquatro.evolve_forward_Euler(1e-4, 100, false, 0, 0, true, "ninetwo_EU_te4");
    ntcinco.evolve_RK4(1e-3, 100, true, 0, 0, true, "ninetwo_RK4_te3_inter");
    */
    /*
    Particle ntto(charge, mass, {1,0,1}, {0,1,0});
    Particle nttt(charge, mass, {-1,0,-1}, {1,0,1});

    PenningTrap nttuno(B0, V0, d);
    PenningTrap ntttres(B0, V0, d);

    nttuno.add_particle(ntto);
    nttuno.add_particle(nttt);
    ntttres.add_particle(ntto);
    ntttres.add_particle(nttt);

    nttuno.evolve_RK4(1e-3, 100, false, 0, 0, true, "ninethree_RK4_te3");
    ntttres.evolve_RK4(1e-3,100,true,0,0,true,"ninethree_RK4_te3_inter");
    */

    /*
    double d_10 = 0.05*cm;
    double V0_10 = 0.0025*volt;
    arma::vec f = {0.1,0.4,0.7};
    
    PenningTrap ten(B0, V0_10, d_10);
    ten.add_n_random_particles(100, charge, mass);
    ten.evolve_RK4(1e-3, 500, f(0))
    */
    /*
    PenningTrap random_test(B0, V0_10, d_10);
    random_test.add_n_random_particles(5, charge, mass);
    random_test.evolve_RK4(1e-3, 100, true, 0,0,true, "RK4_random_test_seed_1");
    std::cout << "Particles left in the penningtrap: " << random_test.particles_inside_trap_count() << "\n";
    */

    return 0;
}