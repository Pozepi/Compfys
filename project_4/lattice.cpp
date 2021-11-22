#include "lattice.hpp"

Lattice::Lattice(int L, double T, bool ordered)
{
    /*
    CONSTRUCTOR

    Creates an instance of the Lattice class.
    Args:
        L       (double)    :   The size of the lattice (LxL).
        T       (double)    :   Temperature of the system.
        ordered (bool)      :   If lattice elements should be ordered or random.
    */
    L_ = L;
    N_ = L_*L_;
    T_ = T;
    
    arma::arma_rng::set_seed_random();
    lattice = arma::randi(L_, L_, arma::distr_param(0,1))*2 - 1; 

    E = Lattice::Total_energy();
    M = Lattice::Total_magnetization();

    double E0 =  8; double expE0 = exp(-E0/T_);
    double E1 =  4; double expE1 = exp(-E1/T_);
    double E2 = -0; double expE2 = exp(-E2/T_);
    double E3 = -4; double expE3 = exp(-E3/T_);
    double E4 = -8; double expE4 = exp(-E4/T_);

    my_map = {
    { E0, expE0},
    { E1, expE1}, 
    { E2, expE2},
    { E3, expE3},
    { E4, expE4}
    }; 
}

arma::mat Lattice::Create_lattice()
{
    /*
    Creates an empty matrix of size (L_ x L_).
    Return:
        lat (arma::mat) :   empty matrix representing lattice.
    */
    arma::mat lat(L_,L_);
    return lat;
}

void Lattice::Fill_lattice(bool ordered)
{   
    /*
    Fills the lattice with values.
    Args:
        ordered (bool)  :   If true, sets all elements in lattice to 1. If false, sets all elements randomly to -1 or 1.
    */
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
    /*
    Expands the lattice with an edge representing the periodic boundary conditions.
    Args:
        lat     (arma::mat) :   lattice of size (L x L) which will be expanded.
    Returns:
        padded  (arma::mat) :   expanded lattice of size (L+2 x L+2).
    */
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
    /*
    Removes the edge of the extended matrix.
    Args:
        pad         (arma::mat) :   matrix (L+2 x L+2) which will have its edge removed.
    Returns:
        unpadded    (arma::mat) :   matrix of size (L x L) where the edge has been removed.
    */
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
    /*
    Calculates the total magnetization of the lattice.
    Returns:
        sum     (double)    :   total magnetization of the lattice.
    */
    double sum = arma::accu(lattice);
    return sum;
}

double Lattice::Total_energy()
{
    /*
    Calculates the total energy of the lattice.
    Returns:
        -sum    (double)    :   total energy of the lattice.
    */
    double sum = 0;

    for(int i=0; i<L_; i++)
    {
        for(int j=0; j<L_; j++)
        {
            sum += lattice(i,j)*(lattice(periodic(i,-1),j) + lattice(i,periodic(j,-1)));
        }
    }
    return -sum;
}
//////////////////////////////////////////
arma::mat Lattice::Replace_pad(arma::mat pad)
{
    /*
    Replaces the elements at the edge of an extended matrix with inner values.
    Args:
        pad     (arma::mat)     :   The matrix whose edge shall be replaced.
    Returns:
        new_pad (arma::mat)     :   A matrix whose edges are correct relative to the inner elements.
    */
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
    /*
    Calculates the energy per spin using the total_energy() function.
    Returns:
        eps     (double)    :   Energy per spin
    */
    double eps = Total_energy()/N_;
    return eps;
}

double Lattice::energy_per_spin_expectation(arma::vec average)
{
    /*
    Calculates the expectation value for the energy per spin.
    Args:
        average         (arma::vec)     :   A vector containing the energy per spin at (6) and number of cycles at (5).
    Returns:
        first_moment    (double)        :   The expectation value for the energy per spin. 
    */  
    double first_moment = average(6)/average(5);
    return first_moment;
}

double Lattice::magnetization_per_spin()
{
    /*
    Calculates magnetization per spin using the total_magnetization() function.
    Returns:
        m   (double)    :   magnetization per spin.
    */
    double m = Total_magnetization()/N_;
    return m;
}

double Lattice::magnetization_per_spin_expectation(arma::vec average)
{
    /*
    Calculates the expectation value for the magnetization per spin.
    Args:
        average         (arma::vec) :   A vector containing the magnetization per spin at (7) and number of cycles at (5).
    Returns:
        first_moment    (double)    :   The expecation value for magnetization per spin. 
    */
    double first_moment = average(7)/average(5);
    return first_moment;
}

double Lattice::specific_heat_capacity(arma::vec average)
{
    /*
    Calculates the specific heat capacity.
    Args:
        average (arma::vec) :   A vector containing the energy at (0), energy^2 at (1) and number of cycles at (5).
    Returns:
        Cv      (double)    :   Specific heat capacity.
    */
    double first_moment = average(0)/average(5);
    double second_moment = average(1)/average(5);
    double Cv = (second_moment - first_moment*first_moment)/(N_*T_*T_);
    return Cv;
}

double Lattice::susceptibility(arma::vec average)
{
    /*
    Calculates the susceptibility.
    Args:
        average (arma::vec) :   A vector containing the absolute value of the total magnetization at (4), total magnetization^2 at (3) and number of cycles at (5).
    Returns:
        chi     (double)    :   Susceptibility.
    */
    double first_moment = average(4)/average(5);
    double second_moment = average(3)/average(5);

    double chi = (second_moment-first_moment*first_moment)/(N_*T_);
    return chi;
}

int Lattice::periodic(int i, int add) 
{
    /*
    Finds the index of the next neighbour of element i in direction add.
    Args:
        i                   (int)   :   index number
        add                 (int)   :   direction
    Returns:
        (i+L_+add) % (L_)   (int)   :   Index of the neighbour.
    */
    return (i+L_+add) % (L_);
}

bool Lattice::test_flip(int i,int j)
{
    /*
    Tests if we should accept the spin flip on index i,j.
    Args:
        i       (int)   :   first index of lattice.
        j       (int)   :   second index of lattice.
    Returns:
        true    (bool)  :   if the new state was accepted.
        false   (bool)  :   if the new state was rejected.
    */
    dE = 2*lattice(i,j)*
    (
     lattice(i,periodic(j,-1))+
     lattice(i,periodic(j,+1))+ 
     lattice(periodic(i,-1),j)+
     lattice(periodic(i,+1),j)
    );

    r = ((double) rand() / (RAND_MAX));
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
    /*
    Runs the metropolis algorithm once.
    Args:
        average (arma::vec&)    :   A vector of results to be filled out.
    */
    for(int k = 0; k < N_; k++)
    {
        int i = arma::randi<int>( arma::distr_param(0,L_-1) );
        int j = arma::randi<int>( arma::distr_param(0,L_-1) );
        if (test_flip(i,j)) // then it passed and we update
        {
            lattice(i,j) *= -1; // flip spin
            E += dE;
            M += 2*lattice(i,j);
        }
    }

    average(0) += E;
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
    /*
    Runs the one_cycle_MCMC() function 50 000 times first to stabilize the system without sampling. 
    Then runs the function again, but sampling second time. 
    Args:
        cycles              (int)           :   Number of times the one_cycle_MCMC() function shall be ran.
        eps_list            (arma::vec&)    :   A vector to be filled out with energy per spin values.
        m_list              (arma::vec&)    :   A vector to be filled out with magnetization per spin values.
        sample_eps_lattice  (bool)          :   If true, fills eps_list and m_list with values directly calculated from the lattice at each cycle.
                                                If false, fills eps_list and m_list with values calculated using the total sum of energy and magnetization per cycle.
    Returns:
        average             (arma::vec)     :   A vector containing results of the samples.  
    */
    int burn_in = 5e4;
    arma::vec average(8);
    average.zeros();
    
    for(int i=0; i<burn_in; i++)
    {
        std::srand((unsigned)time(NULL)+i);
        one_cycle_MCMC(average);
    }
    average.zeros();
    if(sample_eps_lattice==true)
    {
        for(int i = 0; i < cycles; i++)
        {
            std::srand((unsigned)time(NULL)+i);
            one_cycle_MCMC(average);
            eps_list(i) = energy_per_spin();
            m_list(i) = magnetization_per_spin();
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

