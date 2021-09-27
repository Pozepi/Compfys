# main.cpp

main.cpp is a script for solving matrix eigenvalues and eigenvectors analytically and using the Jacobi rotation algorithm. This script requires the functions.cpp and functions.hpp scripts in the same directory in order to run. 

## Compiling and run

Compile and run the script by typing
```
make all
```
in the terminal. This will compile the main.cpp and functions.cpp scripts, create a program named test.out and run the program. If the test.out program already exists in directory, 

# plot_func.py

plot_func.py is a python script for plotting number of iterations as a function of the matrix column size, in addition to plotting a simple model of the expected number of iterations. 

## Installation
Install numpy with 
```
pip install numpy
```

Install matplotlib.pyplot with 
```
pip install matplotlib
```

## Usage 

The script can by run by typing 
```
python plot_func.py
```
In the terminal. 

# plot_last.py

plot_last.py is a python script for plotting the three eigenvectors with the smallest cooresponding eigenvalues, and the analytical values of the eigenvectors. 

## Installation
Install numpy with 
```
pip install numpy
```

Install matplotlib.pyplot with 
```
pip install matplotlib
```

## Usage 
The script can by run by typing
```
python plot_last.py
```
in the terminal window.\\
Files required for this script are:
Markup :
* linspace\_N\_10.txt
* linspace\_N\_100.txt
* N\_10.txt
* N\_100.txt
After plotting the values, the script will also save the figures to files called N\_10.png and N\_100.png.

# functions.cpp

functions.cpp is a c++ script containing functions required to run the main.cpp script.

## Header files required to run the function

Markup :
* armadillo
* iostream
* cmath

## Usage

The functions found in functions.cpp are as follows:\\
* arma::mat tridiagonal_matrix(double a, double d, double N)\\
The tridiagonal_matrix() function creates and returns an NxN matrix with a-values on the subdiagonal and superdiagonal, and d-values on the main diagonal.\\
* arma::vec eigenvalues(double a, double d, int n)\\
The eigenvalues() functions finds the eigenvalues of a nxn tridiagonal matrix with a-values on the subdiagonal and superdiagonal, and d-values on the main diagonal, analytically, and returns these values.\\
* arma::mat eigenvectors(int n)\\
The eigenvectors() function finds the eigenvectors of a nxn tridiagonal matrix analytically, and returns these vectors.\\
* double largest\_off\_element(arma::mat A, int&k, int& l, int n)\\
largest\_off\_element() finds the value of the largest off diagonal element in a matrix, and returns this values. In addition, the function will also save the index of this largest element into variables, k and l. \\
* void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l, int n)\\
jacobi\_rotate() is a function to perform the Jacobi rotation algorithm on matrix A, and store the eigenvectors into a matrix R. \\
* void jacobi_eigensolver(arma::mat& A, arma::vec& eigenvalues, arma::mat& eigenvectors, const in maxiter, int& counter bool& converged)\\
jacobi_eigensolver() is a function which runs the jacobi\_rotate() function untill the input matrix only consists of values along the main diagonal, and zeros everywhere else, or at least values close to zeros. It will then also save the eigenvalues of the resulting matrix.\\
* void write\_eig\_to\_file(arma::vec V, arma::mat R, int n, std::string filename)\\
write\_eig\_to\_file() is a function for writing the the eigenvalues and cooresponding eigenvectors to a file called [filename].txt. The function will write down the list of eigenvectors and put the cooresponding eigenvalue at the top of this list. \\
* void write\_lin\_to\_file(arma::vec x, int n, std::string filename)\\
write\_lin\_to\_file() is a function for writing down the evenly spaced out x-values used in the script into a file called [filename].txt.
