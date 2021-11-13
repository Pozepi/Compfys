#include <iostream>
#include "lattice.hpp"

int main()
{
    Lattice test(3, 1);
    test.Fill_lattice();
    clock_t t1 = clock();
    for(int i = 0; i < 10; i++)
        test.one_cycle_MCMC(10);
    clock_t t2 = clock();
    double duration_seconds = ((double) (t2-t1))/CLOCKS_PER_SEC;
    std::cout << duration_seconds << '\n';

    return 0;
}



