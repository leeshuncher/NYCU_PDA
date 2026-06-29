CXX      = g++
CXXFLAGS = -O3 -DNDEBUG -std=c++11 -Wall -Wextra
TARGET   = Lab2
SRCDIR   = src
SRCS     = $(SRCDIR)/main.cpp $(SRCDIR)/floorplan.cpp

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SRCS) $(SRCDIR)/floorplan.h
	$(CXX) $(CXXFLAGS) -o $@ $(SRCS)

clean:
	rm -f $(TARGET)
