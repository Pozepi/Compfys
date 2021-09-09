
#ifndef Functions_hpp
#define Functions_hpp
#include <iostream>
#include <math.h>
#include <fstream>
#include <algorithm>
// Include the string library
#include <string>

double * linspace(double, double, double*, int);
double * Poisson_analytical(double, double*, int);
void print_array(double, int);
double * fill_array(double*, int, double);
double * solving_matrix(double, double, double, double, double*, int);
void write_to_file(std::string, double, double, int);
double * Delta(double, double, double*, int); 
double * epsilon(double, double, double*, int); 

#endif