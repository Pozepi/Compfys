all: compile execute

compile: 
	g++ *.cpp -std=c++11 -larmadillo -llapack -lblas -fopenmp -o main.exe

execute:
	./main.exe