# PDA Programming Assignments

Programming assignment repository for the Spring 2026 Physical Design
Automation course by Prof. Mark Po-Hung Lin.

## Repository Layout

| Directory | Assignment | Executable | Notes |
| --- | --- | --- | --- |
| `Partition/` | Hypergraph partitioning | `Lab1` | Three-way balanced partitioner. |
| `FloorPlan/` | Symmetry-aware floorplanning | `Lab2` | B*-tree packing with simulated annealing. |
| `Placement/` | Standard-cell legalization | `Legalizer` | Row/site legalization with density refinement. |
| `RMST/` | Rectilinear MST routing | `RMST` | Sparse Manhattan-octant graph plus Kruskal. |

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

## Usage

```sh
cd Partition
./Lab1 <input_file> <output_file>

cd FloorPlan
./Lab2 <input_file> <output_file>

cd Placement
./Legalizer <alpha> <threshold> <input_file> <output_file>

cd RMST
./RMST <input.dat> <output.dat>
```

## Code Style

Run `clang-format -i` on changed C++ files before committing. The shared style
is defined in `.clang-format`.
