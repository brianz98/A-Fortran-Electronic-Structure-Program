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
- [ ] Loop optimisations - using permutational symmetry
- [ ] Adapting Hartree-Fock with symmetry
- [ ] Ditto for MP2
- [ ] Ditto for CCSD
- [ ] Ditto for CCSD(T)
- [ ] MPI parallelisation
- [ ] GPU adaptation for MP2 and CCSD/(T)
