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

# penning_trap.cpp

penning_trap.cpp is a c++ script containing the PenningTrap class. This class requires the Particle class to run. 

## Header files required to run the script

Markup :
* particle.hpp