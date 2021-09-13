#ifndef Functions_hpp
#define Functions_hpp
#include <iostream>
#include <math.h>
#include <fstream>
#include <algorithm>
// Include the string library
#include <string>

void linspace(double*, double, double, int);
void Poisson_analytical(double*, double*, int);

void print_array(double*, int);
void fill_array(double*, int, double);
void solving_matrix(double*, double*, double*, double*, double*, int);
void solving_special(double*, double*, int);
void write_to_file(std::string, double*, double*, int);
void Delta(double*, double*, double*, int); 
void epsilon(double*, double*, double*, int); 

#endif