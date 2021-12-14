# Compile and run

To compile and run the c++ scripts write
```
make all
```

# main.cpp

main.cpp is a c++ script which runs all the simulations from the report. The script requires the particle.cpp and particle.hpp scripts to run. 

# particle.cpp

particle.cpp is a c++ script containing the necessary functions to run the simulations. This script requires particle.hpp to run. 

## Usage
Create an instance of the Particle class by calling 
```
Particle {instance_name}(double h_, double dt_, double T_, double xc_, double yc_, double sigmax_, double sigmay_, double px_, double py_, int pot_type_, double v0_, std::string name_)
```
where {instance_name} is a self defined name of the instance. The parameters are:
Markup : 
* h_ = spacial resolution
* dt_ = size of timestep
* T_ = Final time
* xc_ = center x coordinate of perturbation
* yc_ = center y coordinate of perturbation
* sigmax_ = width of perturbation in x
* sigmay_ = width of perturbation in y
* px_ = momenta of wave packet in x
* py_ = momenta of wave packet in y
* pot\_type\_ = type of potential (0=no slit, 1=one slit, ...)
* v0_ = strength of potential
* name_ = additional name to distinguish saved files

The functions found in particle.cpp are as follows:

* void Particle::initial_state()

Constructs the initial state of the wavefunction matrix using the initial values given. 

* int Particle::transform_index(int i, int j)

Transforms indexes int i, int j of 2D array to a 1D array

* std::tuple<arma::sp_cx_mat, arma::sp_cx_mat> Particle::construct_AB()

Constructs the A and B matrices as defined in the project

* void Particle::update_system()

Updates the wavefunction matrix by one step

* void Particle::simulate_system()

Simulates the system for T_/dt_ timesteps.

* void Particle::potential(int slits)

Creates the potential matrix with 0-3 slits. 

# plot\_func.py

plot\_func.py is a python script for plotting the results generated in main.cpp.

## Usage

Run the script with 
```
python plot\_func.py
```
You will then be given a couple of choices:
Markup:
* [1]: Plots the probability deviation and creates animations of the waves. 
* [2]: Plots the wave at three different times in a subplot. Also plots only the real and immaginary part of this wave in a similar matter. Additionally, creates animation of this wave. 
* [3]: Plots an interference pattern
* [4]: Plots Quantum tunneling
