CXX = g++
CXXFLAGS = -O2 -Wall
LIBS = -lgmp -lntl

# Targets
all: sophie_germain_primes_aks sophie_germain_primes

sophie_germain_primes_aks: sophie_germain_primes_aks.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LIBS)

sophie_germain_primes: sophie_germain_primes.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LIBS)

clean:
	rm -f sophie_germain_primes_aks sophie_germain_primes
