#!/usr/bin/env python
# coding: utf-8

import psi4
import numpy as np

def lowtri_array_2d(end):
    """
    Returns the lower triangular coordinates
    lowtri_array_2d(3) = [[1,1],[2,1],[2,2],[3,1],[3,2],[3,3]]
    """
    idx = []
    for i in range(end+1):
        for j in range(i):
            idx.append([i, j+1])
    return np.array(idx)

def lowtri_array_4d(eri,end):
    """
    Returns the 4d 'lower triangular' coordinates
    lowtri_array_4d(3) = [[1,1,1,1],[2,1,1,1],[2,2,1,1],[2,2,2,1],[3,2],[3,3]]
    """
    idx = []
    for i in range(end):
        for j in range(i+1):
            for k in range(i+1):
                if k<i:
                    lu = k
                else: 
                    lu = j
                for l in range(lu+1):
                    if abs(eri[i,j,k,l]) > 1e-12:
                        idx.append([i+1,j+1,k+1,l+1,eri[i,j,k,l]])
    return np.array(idx)

bl = 1.1 # Angstroms
ang = 104 # Degrees
psi4.set_memory('2000 MB')
mol = psi4.geometry(f"""
    O 
    H 1 {bl}
    H 1 {bl} 2 {ang}
    units = angstrom
    """)
psi4.set_options({'basis':'sto-3g'})
wfn = psi4.core.Wavefunction(mol, psi4.core.BasisSet.build(mol))
mints = psi4.core.MintsHelper(wfn.basisset())
S = mints.ao_overlap().nph
T = mints.ao_kinetic().nph
V = mints.ao_potential().nph
eri = mints.ao_eri().nph

nbasis = S[0].shape[0]
tril_array = lowtri_array_2d(nbasis)
indices = ((tril_array.T-1)[0],(tril_array.T-1)[1])
s_dat = np.hstack((tril_array,np.array([S[0][indices]]).T))
t_dat = np.hstack((tril_array,np.array([T[0][indices]]).T))
v_dat = np.hstack((tril_array,np.array([V[0][indices]]).T))
eri_dat = lowtri_array_4d(eri[0],7)

np.savetxt('s.dat', s_dat, ['%1d','%1d','%17.15f'], delimiter='\t')
np.savetxt('t.dat', t_dat, ['%1d','%1d','%17.15f'], delimiter='\t')
np.savetxt('v.dat', v_dat, ['%1d','%1d','%17.15f'], delimiter='\t')
np.savetxt('eri.dat', eri_dat, ['%1d','%1d','%1d','%1d','%17.15f'], delimiter='\t')