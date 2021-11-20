all: compile execute

compile: 
	g++-11 *.cpp -std=c++11 -larmadillo -llapack -lblas -fopenmp -o main.exe

execute:
	./main.exe