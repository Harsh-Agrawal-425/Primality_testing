CXX = mpic++
CXXFLAGS = -O2 -Wall
LIBS = -lgmp -lntl

TARGET = parallel_primality
SRC = parallel_primality.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC) $(LIBS)

clean:
	rm -f $(TARGET)

