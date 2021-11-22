# Compile and run

To compile and run the c++ scripts write
```
make all
```

Note; if you're running on a mac, compile with
```
g++-11 *.cpp -std=c++11 -larmadillo -llapack -lblas -Xpreprocessor -fopenmp -o main.exe
```
and then run with 
```
./main.exe
```

# main.cpp

main.cpp is a c++ script which runs the Markov Chain Monte Carlo (MCMC) algorithm described in the report. 

## Header files required to run the script

Markup : 
* armadillo
* iostream
* sys/stat.h
* string
* cstdlib
* time.h
* numeric
* random
* algorithm
* map
* time.h
* omp.h
* chrono
* lattice.hpp

Also requires the lattice.cpp script containing the functions. 

## Functions

main.cpp contains one function, except for the int main() function, which is the void loop\_over\_temp() function.

* void loop\_over\_temp(int L, std::string filename, int cycles)
Here, [int L] is an integer number which symbolizes the matrix size. F.ex. L=2 cooresponds to a 2x2 matrix. [int cycles] is the number of cycles the simulation shall be ran for. This is in addition to 50 000 extra cycles for burn-in. The function will run the MCMC simulation over a range of temperatures, from T0 = 2.1 to T1 = 2.4, with 100 points in between. The temperature loop is parallized. Also saves the results into files called ["Temp"+filename] for the temperature, ["Cv"+filename] for the specific heat capacity, ["chi"+filename] for the susceptibility, ["expect_eps"+filename] for the expectation values of the energy per spin and ["expect_m"+filename] for the expectation values of the magnetization per spin. Also prints the computational time. If you don't want to parallize the temperature loop, simply comment out line 22 in the main.cpp script. 

## Usage

When the script is compiled and ran, you will be given a couple of choices:

Markup :
* [1] : Runs the MCMC algorithm for L = 2 with 10 000 MC cycles. Loops over a temperature range between -1 to 1, with 100 points between. Saves output to files called ["Temp10.txt"] for the temperature, ["Cv10.txt"] for the specific heat capacity, ["chi10.txt"] for the susceptibility, ["eps10.txt"] for the expectation value of the energy per spin and ["m10.txt"] for the expectation value of the magnetization per spin.
* [2] : Runs void loop\_over\_temp() function with arguments L = 20, filename = "_L20.txt" and cycles = 250000.
* [3] : Runs the void loop\_over\_temp() function with arguments L = 40, filename = "_L40.txt" and cycles = 250000.
* [4] : Runs the void loop\_over\_temp() function with arguments L = 60, filename = "_L60.txt" and cycles = 250000.
* [5] : Runs the void loop\_over\_temp() function with arguments L = 80, filename = "_L80.txt" and cycles = 250000. 
* [6] : Runs the void loop\_over\_temp() function with arguments L = 100, filename = "_L100.txt" and cycles = 250000.
* [7] : Runs the void loop\_over\_temp() function with user defined inputs for L and cycles, filename = "_Ln.txt".
* [8] : Runs four simulations with 100000 cycles each to test the burn in time for a L=20 lattice. First runs the simulation for T = 1 and ordered elements, saves energy per spin and magnetization per spin to text files ["eps_burn_in_test_ordered_T1.txt"] and ["m_burn_in_test_ordered_T1.txt"]. Second simulation is for T = 1 with unordered elements. Energy per spin and magnetization per spin is saved to ["eps_burn_in_test_unordered_T1.txt"] and ["m_burn_in_test_unordered_T1.txt"]. Third simulation is for T = 2.4 with ordered elements. Energy per spin and magnetization per spin is saved to ["eps_burn_in_test_ordered_T24.txt"] and ["m_burn_in_test_ordered_T24.txt"]. Last simulation is for T = 2.4 with unordered elements. Energy per spin and magnetization per spin is saved to ["eps_burn_in_test_unordered_T24.txt"] and ["m_burn_in_test_unordered_T24.txt"]. 
* [9] : Runs the simulation for 250000 cycles, T = 1 and T = 2.4, calculates the approximated energy per spin values for both temperatures and saves these to files ["approximate_eps_T1.txt"] for T = 1 and ["approximate_eps_T24.txt] for T = 2.4

# lattice.cpp
lattice.cpp is a c++ script for running the MCMC algorithm on a lattice. 

## Header files required to run the script

Markup :
* armadillo
* iostream
* sys/stat.h
* string
* cstdlib
* time.h
* numeric
* random
* algorithm
* map
* time.h
* omp.h
* chrono
* lattice.hpp

## Usage
Create an instance of the Lattice class by calling 
```
Lattice {instance_name}(int L, double Temp, bool ordered)
```
where {instance_name} is a self defined name of the instance, [int L] is the size of the lattice, [double Temp] is the temperature of the system and [bool ordered] is a true/false statement. If true, all lattice elements will be filled with the values 1, while if false, all lattice elements will be filled randomly with either 1 or -1. 

The functions found in lattice.cpp are as follows:

* arma::mat Lattice::Create_lattice()

Creates a LxL sized matrix representing the lattice, and returns this matrix.

* Void Lattice::Fill\_lattice(bool ordered)

Fills the lattice with random values of 1 or -1 if [ordered = false], or fills all elements in lattice with 1 if [ordered = true].

* arma::mat Lattice::Pad_lattice(arma::mat lat)

Expands [arma::mat lat] with and edge containing the elements of the edge of the input matrix. Returns padded, which is the extended matrix, as described in the report. 

* arma::mat Lattice::unpad(arma::mat pad)

Removes the edges of the input matrix [arma::mat pad] and makes the matrix smaller.

* double Lattice::Total_magnetization()

Calculates the total magnetization of the lattice.

* double Lattice::Total_energy()

Calculates the total energy of the lattice.

* arma::mat Lattice::Replace_pad(arma::mat pad)

Updates the edge of the input matrix [arma::mat pad] with the inner edge values of the matrix.

* double Lattice::energy\_per\_spin()

Calculates the energy per spin of the lattice, using the Total\_energy() function.

* double Lattice::energy\_per\_spin\_expectation(arma::vec average)

Takes in a vector with values [arma::vec average], and calculates the expectation value of energy per spin. 

* double Lattice::magnetization\_per\_spin()

Calculates the magnetization per spin of the lattice, using the Total\_magnetization() function.

* double Lattice::magnetization\_per\_spin\_expectation(arma::vec average)

Takes in a vector with values [arma::vec average], and calculates the expectation value of the absolute value of the magnetization per spin. 

* double Lattice::specific\_heat\_capacity(arma::vec average)

Takes in a vector with values [arma::vec average] and calculates the specific heat capacity of the system.

* double Lattice::susceptibility(arma::vec average)

Takes in a vector with values [arma::vec average] and calculates the susceptibility of the system.

* int Lattice::periodic(int i, int add)

Finds the index given that we have periodic boundary conditions. [int i] is the index of the matrix and [int add] is the index you add to find the neighbouring element. 

* bool Lattice::test\_flip(int i, int j)

Calculates the change in energy of the matrix on index [int i, int j], and checks if r is less than or equal to the probability calculated with the change of energy. 

* void Lattice::one_cycle_MCMC(arma::vec& average)

Runs one Monte Carlo cycle / the metropolis algorithm, and fills out a vector called [average] with results. 

* arma::vec Lattice::full_cycle(int cycles, arma::vec& eps_list, arma::vec& m_list, bool sample\_eps\_lattice)

* Creates a vector called [average] which will be filled with resulting values. First, runs the one\_cycle\_MCMC() function 50 000 times without sampling to get the lattice into a more probable starting state. Then runs the one\_cycle\_MCMC() function equivalent to [int cycles], stated as a function argument. If [bool sample\_eps\_lattice] is set to true, [arma::vec& eps_list] and [arma::vec& m_list] will be filled out with energy per spin and magnetization values per spin calculated from the lattice at the current cycle. If [bool sample\_eps\_lattice] is set to false, [arma::vec& eps_list] and [arma::vec& m_list] will be filled with energy per spin and magnetization per spin values based on the sum of all energies and magnetizations on all previous cycles ran. Returns a vector, average, which contains all results. 

* arma::vec average:
Markup :
* average(0): Sum of all energies for all cycles
* average(1): Sum of all energies squared for all cycles
* average(2): Sum of all magnetizations for all cycles
* average(3): Sum of all magnetizations squared for all cycles
* average(4): Sum of all absolute values of all magnetizations for all cycles
* average(5): Total number of cycles
* average(6): Sum of all energies per spin for all cycles
* avverage(7): Sum all all absolute values of magnetization per spin for all cycles. 

# plot\_func.py
plot\_func.py is a function for plotting the results generated in main.cpp.

## Installation

Install pyarma with
```
pip install pyarma 
```

Install numpy with 
```
pip install numpy
```

Install matplotlib.pyplot with 
```
pip install matplotlib
```

## Usage

Run the script with 
```
python plot\_func.py
```

You will then be given a couple of choices:
Markup : 
* [1]: Plots the burn in time figures
* [2]: Plots the approximated probability function
* [3]: Plots analytical and numerical solution
* [4]: Plots values for L = 20, 40, 60, 80 and 100 and fits a gaussian fit. 



