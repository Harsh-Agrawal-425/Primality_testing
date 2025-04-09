# Compiler and flags
CXX = g++
CXXFLAGS = -O2 -Wall -std=c++17

# External libraries
LDLIBS = -lgmp -lgmpxx -lntl

# Executable name
TARGET = prime_generator

# Source files
SRCS = main.cpp \
       aks_prime.cpp \
       miller_rabin.cpp \
       rsa_utils.cpp \
       common_utils.cpp

# Header files
HEADERS = aks.h aks_prime.h miller_rabin.h rsa_utils.h common_utils.h

# Object files
OBJS = $(SRCS:.cpp=.o)

# Default rule
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean rule
clean:
	rm -f $(OBJS) $(TARGET)
