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
objects := $(patsubst src/%.cpp,src/%.o,$(wildcard src/*.cpp))

# Set complier flags 
CXXFLAG = -static -fopenmp -O3 -std=c++11 -lm -llapacke -llapack -lblas -Wall
all: $(OUTPUTD)
$(OUTPUTD): $(objects)
	$(CXX) $(objects) -o $(OUTPUTD) $(CXXFLAG) -L $(ARMALIB)


clean:
	rm -f src/*.o $(OUTPUTD)

