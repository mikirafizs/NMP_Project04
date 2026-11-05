# Diagonalization of symmetric tridiagonal matrix (quantum harmonic oscillator problem).

The key concept of this project is to solve Poisson equation numerically with the method of successive overrelaxation. To obtain a symmetric tridiagonal matrix we discretize the equation and replace the second derivative with three-term finite difference formula. In project we compare the convergences of solutions for different omega values.

**Author:** Michal Sikora
**Subject:** Numerical Methods in Physics
**Major and Faculty:** Computer Physics, Faculty of Physics and Applied Computer Science, AGH

## Structure of the project:
- 'src/' - source code (main.cpp, sparse.hpp)
- 'data/' - solutions of the Poisson equation
- 'plots/' - graphs and script for plotting the data (PNG)
- 'report/' - report of a project (PDF)

## Requirements:
- g++ with C++17 support (tested with GCC 15)
- LAPACK (header `lapacke.h`, linked with `-llapacke`);
- Ubuntu/Debian: `sudo apt install liblapack-dev`
- Python 3.8+ (tested with 3.14)
- matplotlib: `sudo apt install python3-matplotlib` or `pip install matplotlib`
- numpy: `sudo apt install python3-numpy` or `pip install numpy`

## System
- OS: Ubuntu 26.04.1 LTS, Linux 7.0.0-38-generic (aarch64)
- Virtualization: QEMU virtual machine
- CPU: 4 cores, 1 thread per core, 64-bit ARM (aarch64), vendor ID reported as Apple
- RAM: 3.3 GiB (plus 3.8 GiB swap)
- Compiler: g++ (Ubuntu 15.2.0-16ubuntu1) 15.2.0
- Python: <3.14.4>
- matplotlib: <3.10.7+dfsg1>
- BLAS: OpenBLAS <0.3.32+ds-5>
- LAPACK: liblapack3 <3.12.1>

## Build and run

From the `src/` directory:

    g++ -std=c++17 -O2 -Wall -Wextra main.cpp -o main 


The program writes the norm files to `../data/`, so run it from `src/`
and make sure the `data/` directory exists.

## Plots

From the project root:

    python3 plots/plots.py