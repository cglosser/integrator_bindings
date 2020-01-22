# Compile the C++ code
weights.o: weights.cpp
	g++ -O3 -Wall -Wextra `pkg-config --cflags eigen3` -c weights.cpp

# Compile the fortran integrator modules & generate mod files
predictor_corrector.o predictor_corrector.mod: predictor_corrector.f90
	gfortran -O3 -Wall -Wextra -c predictor_corrector.f90

# Compile the test program -- the link step MUST include the -lstdc++ flag!
a.out: test.f90 predictor_corrector.o weights.o
	gfortran -O3 -lstdc++ $^ 