#include "particle.hpp"

Particle::Particle(int M_, double h_, double dt_, double T_, 
    double xc_, double yc_, double sigmax_, double sigmay_, double px_, double py_, double v0_)
{
    /*
    CONSTRUCTOR
    Creates an instance of the particle class.
    Args: 
        - M_        (int)           :   size of the full system
        - h_        (double)        :
        - dt_       (double)        :   size of the timestep
        - T_        (double)        :   value of the final time
        - xc_       (double)        :   center x coordinate of perturbation
        - yc_       (double)        :   center y coordinate of perturbation
        - sigmax_   (double)        :   width of perturbation in x
        - sigmy_    (double)        :   width of perturbation in y
        - px_       (double)        :   momenta of wave packet in x
        - py_       (double)        :   momenta of wave packet in y
    */
    //
    //M = M_;
    M = 1/h_;
    h = h_;
    dt = dt_;
    T = T_; 
    xc = xc_;
    yc = yc_;
    sigmax = sigmax_; 
    sigmay = sigmay_;
    px = px_;
    py = py_;
    v0 = v0_;

    // construct potential based on input (SLIT/TUNNELING)
    // v0 = v0_;
    V = arma::cx_mat(M, M);

    std::tie(A, B) = (*this).construct_AB();
    u = arma::cx_vec((M-2)*(M-2));
    (*this).initial_state();
}

void Particle::initial_state()
{
    /*
    Constructs the initial state of the wavefunction matrix using a gaussian
    perutbation with y center yc and y width defined by sigmay and x center
    at xc and x width defined by sigmax. px and py is the wave packet momenta.
    */
    
    arma::vec x = arma::linspace(0,1,M-2);
    arma::vec y = arma::linspace(0,1,M-2);
    double base_x = 2.0*sigmax*sigmax;
    double base_y = 2.0*sigmay*sigmay;
    const std::complex<double> i_imag(0, 1);

    std::complex<double> element;
    std::complex<double> sum = 0;


    for (int i=1; i<M-1; i++)
    {
        for (int j=1; j<M-1; j++)
        {
            double xxc = x(i-1) - xc;
            double yyc = y(j-1) - yc;
            element = std::exp(-std::pow(xxc, 2)/base_x - std::pow(yyc,2)/base_y + i_imag*px*(xxc) + i_imag*py*(yyc));
            sum += std::conj(element)*element;
            u(transform_index(i,j)) = element;
            // std::cout << i << ' ' << j << '\n'; 
        }
    }
    // Normalize u
    for (int k = 0; k<(M-2)*(M-2); k++)
        u(k) = u(k)/std::sqrt(sum);
    /*
    sum = 0;
    for (int k = 0; k<(M-2)*(M-2); k++)
        sum += std::conj(u(k))*u(k);

    std::cout << sum << '\n';
    */
}


int Particle::transform_index(int i, int j)
{
    /*
    Transforms indexes int i, int j of 2D array to a 1D array.
    i, j must be 0 < i,j < M-1.
    Args: 
        i (int)             : index along first axis
        j (int)             : index along second axis
    Returns: 
        i + (M - 2)*j (int) : index in linear space given i and j.
    */
   return (i-1) + (M - 2)*(j-1);
}

std::tuple<arma::sp_cx_mat, arma::sp_cx_mat> Particle::construct_AB()
{
    /*
    Constructs A and B complex matrices used in updating the matrix.
    Returns:
        - std::make_tuple(A, B) (tuple<cx_mat, cx_mat>) : Tuple containing complex matrix A
        and B.
    */
    
    int n = M - 2;
    int N = std::pow(n, 2);

    // define constants. 
    std::complex<double> r = 1i*dt/(2*h*h);
    std::complex<double> constant1 = 1. + 4.*r;
    std::complex<double> constant2 = 1. - 4.*r;
    std::complex<double> constant3 = 1i*dt/2.;
    
    

    arma::cx_vec a(N);
    arma::cx_vec b(N);

    std::complex<double> z1, z2;

    for (int i = 1; i < n+1; i++)
    {
        for (int j = 1; j < n+1; j++)
        {
            z1 = constant1 + constant3*V(i-1,j-1);
            z2 = constant2 - constant3*V(i-1,j-1);

            a(transform_index(i, j)) = z1;
            b(transform_index(i, j)) = z2;

        }
    }

    arma::sp_cx_mat A(N,N);
    arma::sp_cx_mat B(N,N);


    // Handle the diagonal
    for (int i = 0; i < N; i++)
    {
        A(i,i) = a(i);

        B(i,i) = b(i);
    }
    // Next to diagonal
    A(1,0) = -r;
    A(0,1) = -r;

    B(1,0) =  r;
    B(0,1) =  r;
    for (int i = 0; i < N-1; i++)
    {
        if ((i+1)%n != 0)
        {
            A(i+1,i)    = -r;
            A(i,i+1)    = -r;

            B(i+1,i)    =  r;
            B(i,i+1)    =  r;
        }

    }
    // remaining elements
    for (int i = 0; i < N - n; i++)
    {
        A(i+n,i) = -r;
        A(i,i+n) = -r;

        B(i+n,i) =  r; 
        B(i,i+n) =  r;
    }

    return std::make_tuple(A, B);
}

void Particle::update_system()
{
    /*
    Updates the wavefunction matrix by one step.
    */
    arma::cx_vec b = B * u;
    u = arma::spsolve(A, b);
}

void Particle::potential(int slits)
{  
    double thickx = 0.02;
    double posxy = 0.5;
    double v0 = 10e+5;
    double index_start = std::ceil(posxy/h-thickx/(h*2));
    double cols = std::ceil(thickx*(M));

    for(int i = 0; i < cols; i++)
    {
        for(int j = 0; j < M; j++)
        {
            V(j, index_start-std::floor(cols/2)+i) = v0;
        }
    }   
    switch(slits)
    {
        case 0:
        {
            break;
        }
        case 1:
        {
            break;
        }
        case 2:
        {

            double sep = 0.05;
            double index_starty = std::ceil(posxy/h-sep/(h*2));
            double rows = std::ceil(sep*M);

            for(int i = index_start-std::floor(cols/2); i <= index_start+std::floor(cols/2); i++)
            {
                for(int j = 0; j < rows; j++)
                {
                    V(index_start-rows-std::floor(rows/2)+j, i) = 0;
                    V(index_start+std::ceil(rows/2)+j, i) = 0;
                }
            }
            /*
            for(int i = index_start-std::floor(cols/2); i <= index_start+std::floor(cols/2); i++)
            {
                for(int j = 0; j < M; j++)
                {
                    std::cout << V(j, i) << '\n';
                }
                std::cout << "------------------"<< '\n';
            }
            */

            break;
        }
        case 3:
        {
            break;
        }
    }
}

