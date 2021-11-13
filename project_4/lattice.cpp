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
    double g = 0;
    for(int i = 0; i < L_; i++)
    {
        for(int j=0; j<L_; j++)
        {
            double k = ((std::rand()%2)-0.5)*2;
            lattice(i,j) = g;
            g++;
        }
    }
}
/*
arma::mat Lattice::Pad_lattice(arma::mat lat)
{
    arma::mat padded(L_+2, L_+2);
    for(int i=1; i < L_+1; i++)
    {
        for(int j=1; j < L_+1; j++)
        {
            padded(i,j) = lat(i-1,j-1);
            double *p2 = new double;
            *p2 = padded(L_, j);
            padded(0,j) = *p2;
            double *p3 = new double;
            *p3 = padded(i,L_);
            padded(i,0) = *p3;
            double *p4 = new double;
            *p4 = padded(1, j);
            padded(L_+1,j) = *p4;
            double *p5 = new double;
            *p5 = padded(i,1);
            padded(i,L_+1) = *p5;
        }
    }
    double *p6 = new double;
    *p6 = padded(L_,L_);
    padded(0,0) = *p6;
    double *p7 = new double;
    *p7 = padded(1,L_);
    padded(L_+1,0) = *p7;
    double *p8 = new double;
    *p8 = padded(L_,1);
    padded(0,L_+1) = *p8;
    double *p9 = new double;
    *p9 = padded(1,1);
    padded(L_+1,L_+1) = *p9;
    return padded;
}
*/

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
    //funker ikke for matriser L_<3
    return sum;
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
/*
double Lattice::Boltzman(arma::mat lat)
{
    double Z;
    for (i<)
    double P = exp(-Total_energy()/T)/Z;
    return P;
}
*/
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

void Lattice::one_cycle_MCMC(int n)
{
    arma::mat S = lattice;
    arma::mat pad_s = Pad_lattice(S);
    std::srand(time(NULL));

    double E0 = -8; double expE0 = exp(E0/T_);
    double E1 = -6; double expE1 = exp(E1/T_);
    double E2 = -4; double expE2 = exp(E2/T_);
    double E3 = -2; double expE3 = exp(E3/T_);
    double E4 = -0; double expE4 = exp(E4/T_);

    std::map<double, double> my_map = {
    { E0, expE0},
    { E1, expE1},
    { E2, expE2},
    { E3, expE3},
    { E4, expE4}
    }; 
    
    for(int k = 0; k < n; k++)
    {
        arma::mat S_prime = pad_s;
        int i = (std::rand()%L_)+1;
        int j = (std::rand()%L_)+1;
        S_prime(i,j) = -S_prime(i,j);
        arma::mat S_ = Replace_pad(S_prime);
        double dE = S_prime(i,j)*(S_prime(i-1,j) + S_prime(i+1,j) + S_prime(i,j-1) + S_prime(i,j+1));
        double one = 1;
        double p = std::min(one, my_map[dE]);

        double r = ((double) rand() / (RAND_MAX));

        if(r <= p)
        {
            pad_s = S_prime;
        }
    }

    //calculate some values from the new S
}