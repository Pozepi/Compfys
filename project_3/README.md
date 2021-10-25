# Compile and run

To compile and run the c++ scripts write
```
make all
```
Given that the makefile is included in the directory. This creates and runs a script called test.out

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

* arma::vec external_E_field(arma::vec r, double time, double f, double wv)

Calculates the external E field working at position [r] at time [time], with amplitude [f] and frequency [wv], and returns the E field. If the position [r] is outside of the dimensions of the PenningTrap, this E field is set to 0.

* arma::vec external_B_field(arma::vec r)

Find the external magnetic field. If the particle is outside of the PenningTrap dimensions, this is set to 0.

* arma::vec force_particle(int i, int j)

Calculates the force acting on a particle with index [i] from a particle with index [j].

* arma::vec total_force_external(int i, double time, double f, double wv)

Calculates the force working on a particle [i] at time [time] as a result of the external E field with amplitude [f] and frequency [wv] found in the external\_E\_field() function and as a result of the external B field fouind in the external\_B\_field() function.

* arma::vec total_force_particles(int i)

Uses the force\_particle() function to find the force working on particle [i] from all other particles inside the PenningTrap instance. 

* arma::vec total_force(int i, bool interaction=true, double time=0, double f=0, double wv=0)

Sums all forces working on particle [i] at time [time], with amplitude [f] and frequency [wv]. If [interaction] is set to [false], the forces as a result of particle interactions will be ignored. 

* void evolve_RK4(double dt, double time_stop, bool interaction, double f, double wv, bool makefile, std::string filename)

Integrates the position and velocity of the particle(s) using the Runge Kutta 4 method. [dt] is the time step, [time_stop] is the time at which the simulation will stop. [interaction] turns on and off particle interaction with [true/false]. [f] and [wv] sets the amplitude and frequency used in external\_E\_field() function. If [makefile] is set to [true], the function will create three files with names [filename]\_t.txt, [filename]\_pos.txt and [filename]\_v.txt containing the times, positions and velocitites of each particle at each time step. 

* void evolve_forward_Euler(double dt, double time_stop, bool interaction, double f, double wv, bool makefile, std::string filename)

Integrates the position and velocity of the particle(s) using the Euler-Cromer method. [dt] is the time step, [time_stop] is the time at which the simulation will stop. [interaction] turns on and off particle interaction with [true/false]. [f] and [wv] sets the amplitude and frequency used in external\_E\_field() function. If [makefile] is set to [true], the function will create three files with names [filename]\_t.txt, [filename]\_pos.txt and [filename]\_v.txt containing the times, positions and velocitites of each particle at each time step. 

* double particles_inside_trap_count()

Counts how many particles are inside the dimensions of the PenningTrap instance. 

# main.cpp

main.cpp is a c++ script used to generate results using the particle.cpp and penning_trap.cpp scripts found above. 

## Files required to run the script

Markup :
* particle.hpp
* particle.cpp
* penning_trap.hpp
* penning_trap.cpp

## Usage

This scripts purpose is solely to test different instances of the PenningTrap class found in penning_trap.cpp, for different initial values, amount of particles and so on. 

# plot_penning.py

plot_penning.py is a python script to handle different result files from the PenningTrap and Particle classes, and for plotting results. 

## Installation 

Install matplotlib.pyplot with 
```
pip install matplotlib
```

Install numpy with
```
pip install numpy
```

Install pyarma with
```
pip install pyarma
```

Install glob with
```
pip install glob2
```

## Usage

The script can be run with 
```
python plot_penning.py
```

The functions contained inside the plot_penning.py script are as follows:

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

* relerror(plot=False, save=False)

If [plot] is set to [true], the function will plot the relative error of Runge Kutta 4 and Euler Cromer. If [save] is set to [True], the plot will be saved as a .pdf.

* relerror2()

Calculates the error convergence rate. 

* plot_vel(filename,plot=False,save=False)

If [plot] is set to [True], the function  plots the phase space with values from a file called [filename].txt. Plots each phase space component (x,y,z) in a subplot. If [save] is set to [true], the plots are saved into a .pdf. 

* plot_multiple_particles(filename, animate=False, anisave=False, plot3d=False, saveplot3d=False,
    plot=False, saveplot=False)

Takes in the file [filename], sends it to the find\_values\_single() function, then assigns each particle their correct positions and velocities from the file. If [animate] is set to [True], these particles will be animated in a 3D-plot. If [anisave] is set to [True], and [animate] is set to [True], the animation is saved to a .mp4 file. If [plot3d] is set to [True] particles are plotted in a 3D-plot. If [saveplot3d] is [True], and [plot3d] is set to [True], the 3D particle plot will be saved to a .pdf file. If [plot] is set to [True] the function creates a 2D plot in the x-y direction of the particles, and this plot is saved into a .pdf file if [saveplot] is set to [True]. 

* plot_particle_count(plot=False, save=False)

If [plot] is set to [True], the function plots the number of particles as a function of the frequency, and saves the plot into a .pdf if [save] is set to [True]. 