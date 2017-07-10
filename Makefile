# 1-DDNLS Project Symplectic Solver
#
# Claudia Castro Castro
# @ SMU
# Fall 2015
#
#  g++ version *** 4:4.6.4-1ubuntu5 0

CXX      = g++ 
CXXFLAGS = -O3 -std=c++0x -p -fopenmp
#CXXFLAGS = -O -std=c++11


DEPEND1 = driver_ddnls_nonlinear_pert_ic.cpp matrix.cpp abc4.cpp 

# executable targets
all : driver.exe

driver.exe : $(DEPEND1)
	$(CXX) $(CXXFLAGS) -o $@ $^

# utilities
clean :
	\rm -rf *.txt *.dat *.exe *.o *~

realclean : clean
	\rm -f *.exe *.o *~
