# PDA Programming Assignments

Programming assignment repository for the Spring 2026 Physical Design
Automation course by Prof. Mark Po-Hung Lin.

## Contents

- `Partition/` - hypergraph partitioning implementation.
- `FloorPlan/` - floorplanning with symmetry constraints.
- `Placement/` - standard-cell legalization.
- `RMST/` - rectilinear minimum spanning tree routing.

Each assignment directory contains its own `Makefile`, source files, problem
statement, and sample data where applicable.

## Build

Build an assignment from its directory:

```sh
cd <assignment>
make
```

Clean generated binaries with:

```sh
make clean
```
