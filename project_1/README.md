
# main.cpp

main.cpp is a c++ script for solving the Poisson equation analytically and numerically, and finding the error and relative error. This script requires the script functions.cpp and the header file functions.hpp to exists in the same directory. 

## Compiling
If a main.out file does not exist in your directory, compile the script with this command in your terminal:
```
g++ main.cpp functions.cpp -o main.out
```

## Usage
Run the main.out script with
```
./main.out
```
in the terminal window. You will then have four options to chose from:
Option [1]: Solves the analytical solution of the Poisson equation, and writes the output values to a file called values1.txt.
Option [2]: Solves the numerical solution of the Poisson equation for n=10, 100, 1000, 10000 and writes output values to files called u\_values\_N\_{0,1,2,3}.txt.
Option [3]: Find the errors and relative errors between the analytical and numerical solutions for n=10, 100, 1000, 10000, and writes output values to error\_values\_N\_{0,1,2,3}.txt and rel\_error\_values\_N\_{0,1,2,3}.txt files. 
Option [4]: Finds the maximum relative error for n=10, 100, 1000, 10000, 100000, and writes output values to max\_eps.txt file. 
Option [5]: Times the general and special algorithm for n=10, 100, 1000, 10000, 100000, and writes the time taken for these algorithms to run into time\_general.txt and time\_special.txt.

# plot_func.py

plot\_func.py is a python script for plotting output values from the c++ scripts in this folder. The .txt output files from the c++ scripts need to exist in the same folder as plot\_func.py

## Installation
Install matplotlib with
```
pip install matplotlib
```

Install numpy with 
```
pip install numpy
```

## Usage
To run this program, write
```
python plot_func.py
```
in the terminal window. You will then have four options:
Option [1]: Plot the analytical Poisson solution, with values contained in the values.txt file. 
Option [2]: Plot the numeric Poisson values for n=10, 100, 1000, 10000, found in the four u\_values\_N\_{0,1,2,3}.txt files against the analytical values from values.txt.
Option [3]: Will plot the error and relative errors for n=10, 100, 1000, 10000 contained in the four error\_values\_N\_{0,1,2,3}.txt files and in the four rel\_error\_values\_N\_{0,1,2,3}.txt files. 
Option [4]: Plot the max relative error values for n=10, 100, 1000, 10000, 100000 contained in the max\_eps.txt file. 
Option [5]: Plots the times taken for the general and special algorithm with values found in time\_general.txt and time\_special.txt files. 

# functions.cpp

functions.cpp is a c++ script containing functions required to run the main.cpp script.

## Header files required to run the functions
Markup :
* iostream
* math.h
* fstream
* algorithm
* string

## Usage
The functions found in functions.cpp are as follows:
Markup :
* void linspace(double dummy[], double start, double end, int n)
The linspace() function creates an array with [n] evenly spaced out values from [start] to [end]. An empty array, [dummy], with [n] length must be provided to be fille out. 
Markup :
* void Poisson\_analytical(double dummy[], double x[], int n)
The Poisson\_analytical() function solves the analytical Poisson equation. An empty [dummy] array with [n] length must be provided, as well as an array, [x], with [n] length containing values to base the calculations on. 
Markup :
* void print\_array(double array[], int n)
print\_array() prints a provided array [array] with [n] length.
Markup :
* void fill\_array(double dummy[], int n, double element)
fill\_array() fills the [dummy] array of size [n] with an element provided in the [element] argument.
Markup :
* void solving\_matrix(double dummy[], double x[], double a\_array[], double b\_array[], double c\_array[], int n)
solving\_matrix() solves a matrix A to be used in the numerical solution of the Poisson equation. An empty [dummy] array with length [n] must be provided, as well as an array [x] containing x-values and [a\_array], [b\_array] and [c\_array] containing the values for the subdiagonal, main diagonal and superdiagonal of matrix A. 
Markup :
* void solving\_special(double dummy[], double x[], int n)
solving\_special() solves the numerical solution of the Poisson equation with a specialized algorithm, which is faster than the general algorithm found in the solving\_matrix() function. Empty [dummy] array with length [n] is provided and filled out based on x-values found in [x].
Markup :
* void write\_to\_file(std::string filename, double array[], double lin[], int n)
write\_to\_file() is a function for writing two arrays, [array] and [lin], with length [n] into a file with filename [filename].txt. 
Markup :
* void Delta(double dummy[], double u[], double v[], int n)
Delta() is a function for calculating the absolute error between arrays [u] and [v] with length [n]. An empty array, [dummy], with length [n] must be provided to be filled out. 
Markup :
* void epsilon(double dummy[], double u[], double v[], int n)
epsilon() is a function for calculating the relative error between arrays [u] and [v] with length [n]. An empty array, [dummy], with length [n] must be provided to be filled out.

