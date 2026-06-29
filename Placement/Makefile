CXX := g++
CXXFLAGS := -O3 -std=c++17 -Wall -Wextra -pedantic
TARGET := Legalizer
SRC := src/main.cpp src/Legalizer.cpp

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $@ $(SRC)

clean:
	rm -f $(TARGET)
