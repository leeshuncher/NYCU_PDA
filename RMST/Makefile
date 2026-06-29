CXX ?= g++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -pedantic -DNDEBUG

TARGET := RMST
SRC := src/main.cpp

.PHONY: all clean test

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $@ $<

test: $(TARGET)
	python3 tests/test_rmst.py

clean:
	rm -f $(TARGET) .tmp_*.out
