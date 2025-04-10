CXX = g++
CXXFLAGS = -O2 -Wall
LIBS = -lgmp -lntl

TARGET = data_generation
SRC = data_generation.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC) $(LIBS)

clean:
	rm -f $(TARGET)
