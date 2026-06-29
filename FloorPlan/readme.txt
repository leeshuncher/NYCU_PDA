FloorPlan

This project is a C++ floorplanning program for placing rectangular blocks
with symmetry constraints.

Files
-----
- src/main.cpp: program entry point
- src/floorplan.cpp, src/floorplan.h: floorplanning logic
- Makefile: build script
- *.in: sample input files

Build
-----
Run:

    make

This builds the executable:

    ./Lab2

To remove the executable:

    make clean

Usage
-----
Run the program with an input file and an output file:

    ./Lab2 <input_file> <output_file>

Example:

    ./Lab2 block_30.in block_30.out

The program reads the floorplan problem from block_30.in and writes the block
placements to block_30.out.

Input Format
------------
The input file starts with the number of blocks:

    NumBlocks: 30

Each following block line contains:

    <block_name> <width> <height>

Symmetry constraints are written in groups:

    Symmetry Group
    b1 b2
    b3 b4
    b5

Two names on a line mean the two blocks are a symmetry pair.
One name on a line means the block is self-symmetric in that group.

Output Format
-------------
Each output line contains:

    <block_name> <x> <y> <rotation>

rotation is:

    0 = original width and height
    1 = rotated by 90 degrees

Example Workflow
----------------
1. Build the project:

       make

2. Run a sample case:

       ./Lab2 block_30.in block_30.out
