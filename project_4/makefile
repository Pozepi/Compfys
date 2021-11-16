all: compile execute

compile: 
	g++ *.cpp -std=c++11 -larmadillo -llapack -lblas -fopenmp -o test.out

execute:
	./test.out