#!/usr/bin/env python
# coding: utf-8

import psi4
import numpy as np

bl = 2.0787 # Angstroms
ang = 104.45 # Degrees
psi4.set_memory('2000 MB')
mol = psi4.geometry(f"""
    O 0 -0.143225816552 0
    H 1.638036840407 1.136548822547 0
    H -1.638036840407 1.136548822547 0
    units = bohr
    """)
psi4.set_options({'basis':'sto-3g'})
wfn = psi4.core.Wavefunction(mol, psi4.core.BasisSet.build(mol))
mints = psi4.core.MintsHelper(wfn.basisset())
S = mints.ao_overlap().nph
T = mints.ao_kinetic().nph
V = mints.ao_potential().nph
eri = mints.ao_eri().nph
S_so = mints.so_overlap().nph
T_so = mints.so_kinetic().nph
V_so = mints.so_potential().nph

def lowtri_array(start, interval):
    idx = []
    for i in range(interval+1):
        for j in range(i):
            idx.append([i+start, j+1+start])
    return np.array(idx)

def ints2perm(ints):
    # Input is a tuple of arrays
    nrows = 0
    norb = 0
    for i in range(len(ints)):
        nrows += int(ints[i].shape[0]*(ints[i].shape[0]+1)/2)
    perm = [np.array([]),np.array([]),np.array([])]
    for i in range(len(ints)):
        if len(ints[i]) != 0:
            idx = lowtri_array(norb, ints[i].shape[0])
            norb += ints[i].shape[0]
            perm[0] = np.append(perm[0],idx[:,0])
            perm[1] = np.append(perm[1],idx[:,1])
            perm[2] = np.append(perm[2],lower_triangular(ints[i]))
    return np.stack((perm[:])).T

def lower_triangular(mat):
    flat = []
    for i in range(mat.shape[0]):
        for j in range(i+1):
            flat.append(mat[i,j])
    return flat

s_perm = ints2perm(S_so)

np.savetxt('s.dat', s_perm, ['%1d','%1d','%17.15f'], delimiter='\t')