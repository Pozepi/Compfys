#include <iostream>
#include "monte.hpp"



int main()
{
    Lattice test(3, 1);
    test.Fill_lattice();
    test.one_cycle_MCMC(3);


    return 0;
}



