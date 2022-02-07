# A Fortran Electronic Structure Programme (AFESP)
This project is based on the Crawford Group's excellent [C++ Programming Tutorial in Chemistry](https://github.com/CrawfordGroup/ProgrammingProjects), but written in Fortran. 
The end goal of this project will be performing HF, MP2, CCSD, and CCSD(T), as per the original tutorial, but with additional support for multicore processors (modern CPUs, GPUs).

## Progress / to-do
- [x] Geometry read-in
- [x] Integral read-in
- [x] Nuclear repulsion energy calculation
- [x] Hartree-Fock without symmetry
- [x] Pretty printing, time reporting
- [x] MP2 without symmetry
- [X] CCSD without symmetry
- [X] OpenMP parallelisation
- [X] DIIS acceleration for Hartree-Fock SCF
- [X] CCSD(T) without symmetry
- [X] DIIS acceleration for CCSD iterations
- [X] Loop optimisations - cache
- [X] User input file: level of theory, tolerances
- [X] BLAS acceleration of tensor contractions
- [ ] Spin-free CCSD
- [ ] Spin-free CCSD(T)
- [ ] Loop optimisations - using permutational symmetry
- [ ] Adapting Hartree-Fock with symmetry
- [ ] Ditto for MP2
- [ ] Ditto for CCSD
- [ ] Ditto for CCSD(T)
- [ ] MPI parallelisation
- [ ] GPU adaptation for MP2 and CCSD/(T)

## Installation guide
To get started, clone this directory by
```
git clone git@github.com:brianz98/A-Fortran-Electronic-Structure-Programme.git
```
### Dependencies
You need to have:
- CMake (at least 3.12)
- A Fortran compiler (gfortran 9.3.0 have been tested)
- OpenBLAS (see below for instructions)

### Installing OpenBLAS
Download the latest OpenBLAS version from [their website](https://www.openblas.net/), and untar, then
```
make USE_OPENMP=1
make install [PREFIX=/your/directory]
```
where the `PREFIX` can be omitted and OpenBLAS will be installed in `/opt/OpenBLAS` by default, but you'll likely need to prefix `sudo` to the command.

If you install OpenBLAS in directories other than `/opt/OpenBLAS`, then you need to edit the last line of `CMakeLists.txt` to where `libopenblas.a` actually is.
