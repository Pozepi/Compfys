#include "lattice.hpp"

Lattice::Lattice(int L, double T, bool ordered)
{
    L_ = L;
    N_ = L_*L_;
    T_ = T;
    
    arma::arma_rng::set_seed_random();
    lattice = arma::randi(L_, L_, arma::distr_param(0,1))*2 - 1; //Lattice::Create_lattice();
    // std::cout << lattice << '\n';
    
    //Lattice::Fill_lattice(ordered);

    E = Lattice::Total_energy();
    M = Lattice::Total_magnetization();

    double E0 =  8; double expE0 = exp(-E0/T_);
    double E1 =  4; double expE1 = exp(-E1/T_);
    double E2 = -0; double expE2 = exp(-E2/T_);
    double E3 = -4; double expE3 = exp(-E3/T_);
    double E4 = -8; double expE4 = exp(-E4/T_);
    // std::cout << T_ << '\n'; 

    my_map = {
    { E0, expE0},
    { E1, expE1}, 
    { E2, expE2},
    { E3, expE3},
    { E4, expE4}
    }; 
    // std::cout << my_map[8] << '\n';
}

arma::mat Lattice::Create_lattice()
{
    arma::mat lat(L_,L_);
    return lat;
}

void Lattice::Fill_lattice(bool ordered)
{   
    // arma::mat tmp = arma::randi(L_, L_, arma::distr_param(0,1));

    
    std::srand(time(NULL));
    if(ordered==false)
    {
        for(int i = 0; i < L_; i++)
        {
            for(int j=0; j<L_; j++)
            {
                double k = ((std::rand()%2)-0.5)*2;
                lattice(i,j) = k;
            }
        }
    }
    else if(ordered==true)
    {
        lattice.ones(L_,L_);
    }
    
}
//////////////////////////////////////////
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
//////////////////////////////////////////
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
//////////////////////////////////////////
double Lattice::Total_magnetization()
{
    double sum = arma::accu(lattice);
    return sum;
}

double Lattice::Total_energy()
{
    double sum = 0;

    for(int i=0; i<L_; i++)
    {
        for(int j=0; j<L_; j++)
        {
            sum += lattice(i,j)*(lattice(periodic(i,L_,-1),j) + lattice(i,periodic(j,L_,-1)));
        }
    }
    // std::cout << lattice << '\n';
    // std::cout << "E = " << -sum << '\n';
    return -sum;
}
//////////////////////////////////////////
arma::mat Lattice::Replace_pad(arma::mat pad)
{
    arma::mat new_pad(L_+2, L_+2);
    for(int i=1; i<L_+1; i++)
    {
        for(int j=1; j<L_+1; j++)
        {
            new_pad(i,j) = pad(i,j);
            new_pad(0,j) = pad(L_, j);
            new_pad(i,0) = pad(i, L_);
            new_pad(L_+1,j) = pad(1, j);
            new_pad(i,L_+1) = pad(i, 1);
        }
    }
    new_pad(0,0) = pad(L_,L_);
    new_pad(L_+1,0) = pad(1,L_);
    new_pad(0,L_+1) = pad(L_,1);
    new_pad(L_+1,L_+1) = pad(1,1);
    return new_pad;
}

double Lattice::energy_per_spin()
{
    double eps = Total_energy()/N_;
    return eps;
}

double Lattice::energy_per_spin_expectation(arma::vec average)
{
    double first_moment = average(6)/average(5);
    return first_moment;
}

double Lattice::magnetization_per_spin()
{
    double m = Total_magnetization()/N_;
    return m;
}

double Lattice::magnetization_per_spin_expectation(arma::vec average)
{
    double first_moment = average(7)/average(5);
    return first_moment;
}

double Lattice::specific_heat_capacity(arma::vec average)
{
    double first_moment = average(0)/average(5);
    double second_moment = average(1)/average(5);
    // std::cout << first_moment << '\n';
    // std::cout << second_moment << '\n';
    double Cv = (second_moment - first_moment*first_moment)/(N_*T_*T_);
    return Cv;
}

double Lattice::susceptibility(arma::vec average)
{
    double first_moment = average(4)/average(5);
    double second_moment = average(3)/average(5);

    double chi = (second_moment-first_moment*first_moment)/(N_*T_);
    return chi;
}

int Lattice::periodic(int i, int limit, int add) 
{
    return (i+limit+add) % (limit);
}

bool Lattice::test_flip(int i,int j)
{
    dE = 2*lattice(i,j)*
    (
     lattice(i,periodic(j,L_,-1))+
     lattice(i,periodic(j,L_,+1))+ 
     lattice(periodic(i,L_,-1),j)+
     lattice(periodic(i,L_,+1),j)
    );

    r = ((double) rand() / (RAND_MAX));
    // std::cout << my_map[8] << '\n';
    p = std::min(1.0, my_map[dE]);

    if (r <= p)
    {
        return true;
    }
    else
        return false;
}

void Lattice::one_cycle_MCMC(arma::vec& average)
{
    for(int k = 0; k < N_; k++)
    {
        // arma::mat S_prime = padded;
        // std::srand(time(NULL)+k); //Randomize seed initialization
        
        //int i = (std::rand()%L_); // random index in lattice
        //int j = (std::rand()%L_); // ...
        int i = arma::randi<int>( arma::distr_param(0,L_-1) );
        int j = arma::randi<int>( arma::distr_param(0,L_-1) );
        //std::cout << i << ' ' << j << '\n';  

        // std::cout << i << " " << j << '\n';
        if (test_flip(i,j)) // then it passed and we update
        {
            lattice(i,j) *= -1; // flip spin
            // std::cout << lattice << '\n';
            E += dE;
            M += 2*lattice(i,j);
        }
    }
    //double E = Lattice::Total_energy(pad_s, true);
    //double M = Lattice::Total_magnetization(pad_s, true);
    //double eps = Lattice::energy_per_spin(pad_s, true);
    //double m = Lattice::magnetization_per_spin(pad_s, true);

    average(0) += E;
    // std::cout << average(0) << '\n';
    average(1) += E*E;
    average(2) += M;
    average(3) += M*M;
    average(4) += std::fabs(M);
    average(5) += 1;
    average(6) += E/N_;
    average(7) += std::fabs(M)/N_;
}

arma::vec Lattice::full_cycle(int cycles, arma::vec& eps_list, arma::vec& m_list, bool sample_eps_lattice)
{
    int burn_in = 5e4;
    arma::vec average(8);
    average.zeros();
    
    /*
    for(int i=0; i<burn_in; i++)
    {
        std::srand((unsigned)time(NULL)+i);
        one_cycle_MCMC(average);
    }
    */
    if(sample_eps_lattice==true)
    {
        for(int i = 0; i < cycles; i++)
        {
            std::srand((unsigned)time(NULL)+i);
            one_cycle_MCMC(average);
            eps_list(i) = energy_per_spin();
            m_list(i) = magnetization_per_spin();
            std::cout << lattice << std::endl;
            std::cout << energy_per_spin()<< std::endl;
        }
    }
    else if(sample_eps_lattice==false)
    {
        for(int i = 0; i < cycles; i++)
        {
            std::srand((unsigned)time(NULL)+i);
            one_cycle_MCMC(average);
            eps_list(i) = average(0)/(N_);
            m_list(i) = average(4)/(N_);
        }
    }
    return average;
}

