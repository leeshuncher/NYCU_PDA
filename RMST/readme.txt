Build:
  make

Usage:
  ./RMST input.dat output.dat

The program reads points from input.dat, constructs a sparse Manhattan-octant
candidate graph, runs Kruskal's algorithm, and writes the total rectilinear MST
weight as a single signed 64-bit integer to output.dat.

Optional local test:
  python3 tests/test_rmst.py
