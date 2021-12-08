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
	$(CXX) -c src/main.cpp $(CXXFLAG)

src/%.o: src/%.cpp include/%.hpp
	$(CXX) -c %.cpp $(CXXFLAG)



clean:
	rm -f src/*.o $(OUTPUTD)

