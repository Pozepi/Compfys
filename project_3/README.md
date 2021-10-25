# main.cpp






# particle.cpp

particle.cpp is a c++ script containing the Particle class, which creates a Particle element. 

## Header files required to run the script

Markup : 
* particle.hpp
* armadillo

## Usage
Create an instance of the Particle class by calling 
```
Particle [particle_name](charge, mass, position, velocity)
```
Where [charge] is a double and contains the particles charge, [mass] is a double and contains the particles mass, [position] is an arma::vec and contains the particles initital position in x, y and z-directions (x,y,z), and [velocity] is an arma::vec containing the particles initial velocity in x, y and z (x,y,z).

The functions found in particle.cpp are as follows:

* double Particle::charge()

The Particle::charge() function returns the charge of the particle.

* double Particle::mass()

The Particle::mass() function returns the mass of the particle.

* arma::vec Particle::position()

The Particle::position() function returns the position of the particle in (x,y,z) coordinates.

* arma::vec Particle::velocity()

The Particle::velocity() function returns the velocity of the particle in (x,y,z) coordinates. 

# penning_trap.cpp

penning_trap.cpp is a c++ script containing the PenningTrap class. This class requires the Particle class to run. 

## Header files required to run the script

Markup :
* particle.hpp
* penning_trap.hpp
* armadillo
* isostream
* sys/stat.h
* string
* cstdlib

## Usage

Create an instance of the PenningTrap class by calling 
```
PenningTrap [instance name](magnetic_field, potential, dimension)
```
Where [magnetic_field] is the magnetic field B, [potential] is the electric potential V and [dimension] is the dimension of the PenningTrap.

The functions found in penning_trap.cpp are as follows:

* void add_particle(Particle particle_in)

Adds a Particle object [particle_in] from the Particle class into the PenningTrap instance. 

* void add_n_random_particles(int n, double charge, double mass)

Adds n Particles into the PenningTrap instance, with specified charge and mass. Gives each particle a random initial position and velocity.

* double particle_count()

Counts the number of particles inside the PenningTrap instance, and returns the number.

* bool particle_outside_trap_check(arma::vec r)

Checks if a particle at position [r] is still inside the [dimension] of the PenningTrap instance. 

* arma::vec external_E_field(arma::vec r, double time=0, double f=0, double wv=0)

Calculates the external E field in the PenningTrap and returns the E field. If the particle is outside of the dimensions of the PenningTrap, this E field is set to 0. 

* arma::vec external_B_field(arma::vec r)

Find the external magnetic field. If the particle is outside of the PenningTrap dimensions, this is set to 0.

* arma::vec force_particle(int i, int j)

Calculates the 


void add_particle(Particle particle_in);
    void add_n_random_particles(int n, double charge, double mass);
    double particle_count();
    bool particle_outside_trap_check(arma::vec r);
    arma::vec external_E_field(arma::vec r, double time=0, double f=0, double wv=0);
    arma::vec external_B_field(arma::vec r);
    arma::vec force_particle(int i, int j);
    arma::vec total_force_external(int i, double time=0, double f=0, double wv=0);
    arma::vec total_force_particles(int i); 
    arma::vec total_force(int i, bool interaction=true, double time=0, double f=0, double wv=0);
    void evolve_RK4(double dt, double time_stop, bool interaction=true, double f=0, double wv=0, bool makefile=false, std::string filename="RK4");
    void evolve_forward_Euler(double dt, double time_stop, bool interaction=true, double f=0, double wv=0, bool makefile=false, std::string filename="forward_euler");
    double particles_inside_trap_count();

# plot_penning.py



## Usage

* x(t, v0, x0)

Where [t] is the time, [v0] is the initial velocity in the y-direction and [x0] is the initial position in x-direction. The function calculates the x-position of the particle analytically through time. 

* y(t, v0, x0)

Where [t] is the time, [v0] is the initial velocity in the y-direction and [x0] is the initial position in x-direction. The function calculates the y-position of the particle analytically through time. 

* z(t, z0)

Where [t] is the time, [z0] is the initial position in the z-direction. The function calculates the z-position of the particle analytically through time. 

* find_values_single(filename)

Loads time, position and velocity values from three files called [filename]+_t.txt, [filename]+_pos.txt and [filename]+_v.txt, and returns these values into three arrays. 

* plot_compare_analytical(filename1, plot=False)

Calls the find\_values\_single(filename) function with the [filename1] argument. If [plot] is set to [True], the function will plot the xy-plane into one plot and zt-plane into another. Also plots the analytical solution in the same plots. 

* relerror(plot=False)

