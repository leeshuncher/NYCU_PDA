UNAME_S := $(shell uname -s)

ifneq ($(filter default undefined,$(origin CXX)),)
ifeq ($(UNAME_S),Darwin)
ifneq ($(shell command -v g++-15 2>/dev/null),)
CXX := g++-15
else ifneq ($(shell command -v g++-14 2>/dev/null),)
CXX := g++-14
else
CXX := g++
endif
else
CXX := g++
endif
endif

CXXFLAGS ?= -std=c++20 -O3 -DNDEBUG -flto -march=native -Wall -Wextra -pedantic

SRCS := $(wildcard src/*.cpp)

.PHONY: all clean

all: Lab1

Lab1: $(SRCS)
	$(CXX) $(CXXFLAGS) $(SRCS) -o $@

clean:
	rm -f Lab1
