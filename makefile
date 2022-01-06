# CC=mpicc
CC = g++


sysname=$(uname -n)
UNAME := $(shell uname)

# ifeq ($(UNAME), Darwin)
# 	OPTIONS = -O3 -I/usr/local/include -L/usr/local/lib -lfftw3
# endif
# ifeq ($(UNAME), Linux)
# 	OPTIONS = -O3 -I/opt/fftw/3.3.4/intel/mvapich2_ib/include -L/opt/fftw/3.3.4/intel/mvapich2_ib/lib -lfftw3
# endif

# # for cygwin
# OPTIONS = -O3 -I/usr/include -L/usr/lib -lfftw3

# for Mac OS
OPTIONS = -O3 -lm -std=c++11 -I/usr/local/opt/libomp/include -L/usr/local/opt/libomp/lib -lomp

EXE = mhd.exe 
MAIN = mhd.cpp
FILE01 = initialize.cpp usrinitialize.cpp output.cpp support.cpp \
	ios_support.cpp calcDerivs.cpp timeadvance.cpp mhdrhs.cpp filter.cpp

OBJMAIN = ${MAIN:.cpp=.o}
OBJ01 = ${FILE01:.cpp=.o}
OBJ   = $(OBJMAIN) $(OBJ01)


$(EXE): $(OBJ)
	$(CC) -o $(EXE) $(OBJ) $(OPTIONS) 
$(OBJMAIN):
	$(CC) -c $(MAIN) $(OPTIONS)
$(OBJ01):
	$(CC) -c $(FILE01) $(OPTIONS)

clean:
	rm *.o

cleanData:
	rm *.dat rec log
