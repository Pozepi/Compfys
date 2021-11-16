#include "lattice.hpp"

Lattice::Lattice(int L, double T)
{
    L_ = L;
    N_ = L_*L_;
    T_ = T;
    lattice = Lattice::Create_lattice();
    Lattice::Fill_lattice();
    padded = Lattice::Pad_lattice(lattice);
}

arma::mat Lattice::Create_lattice()
{
    arma::mat lat(L_,L_);
    return lat;
}

void Lattice::Fill_lattice()
{   
    std::srand(time(NULL));
    for(int i = 0; i < L_; i++)
    {
        for(int j=0; j<L_; j++)
        {
            double k = ((std::rand()%2)-0.5)*2;
            lattice(i,j) = k;
        }
    }
    
   //lattice.ones(L_,L_);
}

arma::mat Lattice::Pad_lattice(arma::mat lat)
{
    arma::mat padded(L_+2, L_+2);
    for(int i=1; i < L_+1; i++)
    {
        for(int j=1; j < L_+1; j++)
        {
            padded(i,j) = lat(i-1,j-1);
            padded(0,j) = padded(L_, j);
            padded(i,0) = padded(i, L_);
            padded(L_+1,j) = padded(1, j);
            padded(i,L_+1) = padded(i, 1);
        }
    }
    padded(0,0) = padded(L_,L_);
    padded(L_+1,0) = padded(1,L_);
    padded(0,L_+1) = padded(L_,1);
    padded(L_+1,L_+1) = padded(1,1);
    return padded;
}
arma::mat Lattice::unpad(arma::mat pad)
{
    arma::mat unpadded(L_, L_);
    for(int i = 1; i < L_+1; i++)
    {
        for(int j = 1; j < L_+1; j++)
        {
            unpadded(i-1,j-1) = pad(i, j);
        }
    }
    return unpadded;
}

double Lattice::Total_magnetization(arma::mat lat, bool padded)
{
    arma::mat latt = lat;
    if(padded == true)
    {
        latt = unpad(lat);
    }
    double sum = arma::accu(latt);
    return sum;
}

double Lattice::Total_energy(arma::mat lat, bool padded)
{
    double sum = 0;
    arma::mat pad = lat;

    if(padded == false)
        {
        pad = Pad_lattice(lat);
        }
    for(int i=1; i<L_+1; i++)
    {
        for(int j=1; j<L_+1; j++)
        {
            sum += pad(i,j)*(pad(i-1,j)+pad(i,j-1));
        }
    }
    return -sum;
}

arma::mat Lattice::Replace_pad(arma::mat lat)
{
    arma::mat new_lat(L_+2, L_+2);
    for(int i=1; i<L_+1; i++)
    {
        for(int j=1; j<L_+1; j++)
        {
            new_lat(i,j) = lat(i,j);
            new_lat(0,j) = lat(L_, j);
            new_lat(i,0) = lat(i, L_);
            new_lat(L_+1,j) = lat(1, j);
            new_lat(i,L_+1) = lat(i, 1);
        }
    }
    new_lat(0,0) = lat(L_,L_);
    new_lat(L_+1,0) = lat(1,L_);
    new_lat(0,L_+1) = lat(L_,1);
    new_lat(L_+1,L_+1) = lat(1,1);
    return new_lat;
}

double Lattice::Boltzman()
{
    double Z;

    arma::mat dummy;
    dummy.ones(L_, L_);
    int count = 1;
    // Z = ... manually set Z to just ones
    //std::cout << count << '\n';
    //std::cout << dummy << '\n';

    // for (int k = 0; k<N_*N_; k++)
    while (true)
    {   // for all elements in lattice
        // subtract 2 from the first element
        dummy(0,0) -= 2;
        // the thought is to roll this subtraction
        // over one index if we go under -1

        // basically numpy.where
        for (int j = 0; j<L_; j++)
        { 
            for (int i = 0; i<L_; i++)
            {
                // where and elements is less than -1
                if (dummy(i,j) < -1)
                {
                    // we subtract to the next element
                    if (i<L_-1)
                        dummy(i+1,j) -= 2;
                    // if we are on the border we subtract the
                    // first element in the next row
                    if (i==L_-1)
                        dummy(0,j+1) -= 2;
                    // else if we are in the corner (i,j) = (L,L) we
                    // simply do nothing
                    else if (i==L_-1 | j==L_-1)
                        double null = 0;
                    // finally set current index to 1
                    dummy(i,j) = 1;
                }
            }
        }
        // Z += ...
        count += 1;
        
        //std::cout << count << '\n';
        //std::cout << dummy << '\n';

        if (arma::accu(dummy) == -N_)
            break;
        // PRINT HERE SENPAI UWU
    }
    std::cout << count << '\n';

    /*
    for (i<)
    double P = exp(-Total_energy()/T)/Z;
    return P;
    */
   double ten = 10;
   return ten;
}

double Lattice::energy_per_spin(arma::mat lat, bool padded)
{
    double eps = Total_energy(lat, padded)/N_;
    return eps;
}

double Lattice::magnetization_per_spin(arma::mat lat, bool padded)
{
    double m = Total_magnetization(lat, padded)/N_;
    return m;
}

double Lattice::specific_heat_capacity(arma::vec average)
{
    
    double first_moment = average(0)/average(5);
    std::cout << first_moment << '\n';
    double second_moment = average(1)/average(5);
    std::cout << second_moment << '\n';
    double Cv = (second_moment - first_moment*first_moment)/(N_*T_*T_);
    return Cv;
}

double Lattice::susceptibility(arma::vec average)
{
    double first_moment = average(4)/average(5);
    std::cout << first_moment << '\n';
    double second_moment = average(3)/average(5);
    std::cout << second_moment << '\n';
    double chi = (second_moment-first_moment*first_moment)/(N_*T_);
    return chi;
}

void Lattice::one_cycle_MCMC(arma::vec& average, std::map<double, double> my_map)
{
    arma::mat S = lattice;
    arma::mat pad_s = Pad_lattice(S);

    for(int k = 0; k < N_; k++)
    {
        arma::mat S_prime = pad_s;
        int i = (std::rand()%L_)+1; // random index in lattice
        int j = (std::rand()%L_)+1; // ...
        S_prime(i,j) = -S_prime(i,j); // flip random spin within lattice (pad not affected)
        arma::mat S_ = Replace_pad(S_prime); // update pad
        // calculate the change in energy
        double dE = -S_(i,j)*(S_(i-1,j) + S_(i+1,j) + S_(i,j-1) + S_(i,j+1)) + pad_s(i,j)*(pad_s(i-1,j) + pad_s(i+1,j) + pad_s(i,j-1) + pad_s(i,j+1));
        // double dE1 = Lattice::Total_energy(S_, true) - Lattice::Total_energy(pad_s, true);
        double one = 1;
        double p = std::min(one, my_map[dE]);
        // std::cout << my_map[dE] << " " << dE << '\n'; 

        double r = ((double) rand() / (RAND_MAX));

        if(r <= p)
        {
            pad_s = S_;
            lattice = pad_s;
        }
        // std::cout << pad_s << '\n';
    }
    double E = Lattice::Total_energy(pad_s, true);
    double M = Lattice::Total_magnetization(pad_s, true);
    //std::cout << '\n';
    //std::cout << E << '\n'; 
    //std::cout << M << '\n';

    average(0) += E;
    average(1) += E*E;
    average(2) += M;
    average(3) += M*M;
    average(4) += std::fabs(M);
    average(5) += 1;
    // std::cout << average(0) << '\n';

    //calculate some values from the new S
}
/*
void Lattice::one_cycle_MCMC(arma::vec& average, std::map<double, double> my_map)
{
    arma::mat S = lattice;
    arma::mat pad_s = Pad_lattice(S);

    for(int k = 0; k < N_; k++)
    {
        arma::mat S_prime = pad_s;
        int i = (std::rand()%L_)+1;   // random index in lattice
        int j = (std::rand()%L_)+1;   // ...
        S_prime(i,j) = -S_prime(i,j); // flip random spin within lattice (pad not affected)
        arma::mat S_ = Replace_pad(S_prime); // update pad
        // calculate the change in energy
        double dE = -S_(i,j)*(S_(i-1,j) + S_(i,j-1)) + pad_s(i,j)*(pad_s(i-1,j) + pad_s(i,j-1));
        // std::cout << dE << '\n';
        // double dE1 = Lattice::Total_energy(S_, true) - Lattice::Total_energy(pad_s, true);
        // std::cout << dE1 << '\n';
        double one = 1;
        double p = std::min(one, my_map[dE]);
        // std::cout << my_map[dE] << " " << dE << '\n'; 

        double r = ((double) rand() / (double) RAND_MAX);

        if(r <= p)
        {
            pad_s = S_;
        }
    }
    double E = Lattice::Total_energy(pad_s, true);
    double M = Lattice::Total_magnetization(pad_s, true);
    //std::cout << '\n';
    //std::cout << E << '\n'; 
    //std::cout << pad_s << '\n';
    //std::cout << M << '\n';

    average(0) += E;
    average(1) += E*E;
    average(2) += M;
    average(3) += M*M;
    average(4) += std::fabs(M);
    average(5) += 1;
    // std::cout << average(0) << '\n';

    //calculate some values from the new S
}
*/


arma::vec Lattice::full_cycle(int cycles)
{
    arma::vec average(6);
    average.zeros();

    double E0 =  8; double expE0 = exp(-E0/T_);
    double E1 =  4; double expE1 = exp(-E1/T_);
    double E2 = -0; double expE2 = exp(-E2/T_);
    double E3 = -4; double expE3 = exp(-E3/T_);
    double E4 = -8; double expE4 = exp(-E4/T_);

    std::map<double, double> my_map = {
    { E0, expE0},
    { E1, expE1},
    { E2, expE2},
    { E3, expE3},
    { E4, expE4}
    }; 

    for(int i = 0; i < cycles; i++)
    {
        std::srand((unsigned)time(NULL)+i);
        one_cycle_MCMC(average, my_map);
    }

    return average;
}