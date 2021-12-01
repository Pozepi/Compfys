all: compile execute

compile: 
	g++-11 *.cpp -std=c++11 -fext-numeric-literals -larmadillo -llapack -lblas -fopenmp -o main.exe

execute:
	./main.exe