# # EIGENFLAGS=-DEIGEN_USE_BLAS -DEIGEN_USE_LAPACKE # take out if you don't have blas/lapacke installed
# CPPFLAGS= -DNDEBUG -I${mkEigenInc} # -I. -I/usr/include/eigen3 $(EIGENFLAGS)
# #LDLIBS=-llapacke -lblas# take out if you don't have blas/lapacke installed
# CXXFLAGS= -std=c++17 # -O3 -std=c++20 -Wall -Wextra -Wpedantic -Werror
# LDFLAGS=-O3
# mainFem: mainFem.o
# mainFem.o: mainFem.cpp LinearFiniteElement.hpp gradient.hpp
# clean:
# 	rm -f mainFem mainFem.o
# # Rule to link object files into an executable, ensuring the C++ compiler is used
# %: %.o
# 	$(CXX) $(LDFLAGS) -o $@ $^	$(LDLIBS)

# EIGENFLAGS=-DEIGEN_USE_BLAS -DEIGEN_USE_LAPACKE # take out if you don't have blas/lapacke installed
CPPFLAGS= -DNDEBUG -I${mkEigenInc} # -I. -I/usr/include/eigen3 $(EIGENFLAGS)
#LDLIBS=-llapacke -lblas# take out if you don't have blas/lapacke installed
CXXFLAGS= -std=c++17 # -O3 -std=c++20 -Wall -Wextra -Wpedantic -Werror
LDFLAGS=-O3
prova: prova.o
prova.o: prova.cpp LinearFiniteElement.hpp gradient.hpp
clean:
	rm -f prova prova.o
# Rule to link object files into an executable, ensuring the C++ compiler is used
%: %.o
	$(CXX) $(LDFLAGS) -o $@ $^	$(LDLIBS)