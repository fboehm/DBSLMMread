# -----------------------------------------------------------------
#   Makefile for dbslmmread 
# ---------------------------------------------------------------------

# Set the file type
OUTPUTD = bin/dbslmmread

# Set the path of library
ARMALIB = /net/mulan/home/yasheng/Cpp/arma/lib

# Put C++ complier 
CXX = g++
#CXX = clang++-11

# Set complier flags 
CXXFLAG = -static -fopenmp -O3 -std=c++11 -lm -llapacke -llapack -lblas -Wall
all: $(OUTPUTD)
$(OUTPUTD): src/main.o src/calc_asymptotic_variance.o src/dbslmm.o src/dptr.o src/tobool.o
	$(CXX) src/main.o src/calc_asymptotic_variance.o src/dbslmm.o src/dptr.o src/tobool.o -o $(OUTPUTD) $(CXXFLAG) -L $(ARMALIB)
src/main.o: src/main.cpp
	$(CXX) -c src/main.cpp -o src/main.o

src/calc_asymptotic_variance.o: src/calc_asymptotic_variance.cpp include/calc_asymptotic_variance.h 
	$(CXX) -c src/calc_asymptotic_variance.cpp -o src/calc_asymptotic_variance.o $(CXXFLAG)

src/dbslmm.o: src/dbslmm.cpp
	$(CXX) -c src/dbslmm.cpp -o src/dbslmm.o

src/dtpr.o: src/dtpr.cpp
	$(CXX) -c src/dtpr.cpp -o src/dtpr.o

src/tobool.o: src/tobool.cpp
	$(CXX) -c src/tobool.cpp -o src/tobool.o


clean:
	rm -f src/*.o $(OUTPUTD)

