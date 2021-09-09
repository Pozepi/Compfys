
#ifndef Functions_hpp
#define Functions_hpp
#include <iostream>
#include <math.h>
#include <fstream>

double * linspace(double, double, double*, int);
double * Poisson_analytical(double x[]);
void print_array(double*, int);
double * fill_array(double*, int, double);
double * solving_matrix(double*, double*, double*, double*, double*, int);
void write_to_file(double*, double*, int);
double * Delta(double u[], double v[]); 
double * epsilon(double u[], double v[]); 

#endif