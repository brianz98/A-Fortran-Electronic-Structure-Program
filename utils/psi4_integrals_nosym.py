#!/usr/bin/env python
# coding: utf-8

import psi4
import numpy as np
from pathlib import Path

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

def generate_dat_psi(dirname, mol, wfn):
    Path(f'{dirname}').mkdir(parents=True,exist_ok=True)
    mints = psi4.core.MintsHelper(wfn.basisset())
    geom = np.zeros((mol.natom(),4))
    geom[:,1:] = mol.geometry().nph[0]
    geom[:,0] = [int(mol.charge(_)) for _ in range(mol.natom())]
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
    eri_dat = lowtri_array_4d(eri[0],nbasis)

    np.savetxt(f'{dirname}/geom.dat', geom, ['%1d','%17.15f','%17.15f','%17.15f'], delimiter='\t')
    with open(f'{dirname}/geom.dat', 'r') as f:
        geomcontent = f.readlines()
    geomcontent.insert(0, str(mol.natom())+'\n')
    with open(f'{dirname}/geom.dat', 'w') as f:
        geomcontent = "".join(geomcontent)
        f.write(geomcontent)

    np.savetxt(f'{dirname}/s.dat', s_dat, ['%1d','%1d','%17.15f'], delimiter='\t')
    np.savetxt(f'{dirname}/t.dat', t_dat, ['%1d','%1d','%17.15f'], delimiter='\t')
    np.savetxt(f'{dirname}/v.dat', v_dat, ['%1d','%1d','%17.15f'], delimiter='\t')
    np.savetxt(f'{dirname}/eri.dat', eri_dat, ['%1d','%1d','%1d','%1d','%17.15f'], delimiter='\t')

def generate_water(bl, ang):
    mol = psi4.geometry(f"""
        O 
        H 1 {bl}
        H 1 {bl} 2 {ang}
        noreorient
        symmetry c1
        units = angstrom
        """)
    wfn = psi4.core.Wavefunction(mol, psi4.core.BasisSet.build(mol))
    return mol, wfn

def run_psi4(bl, ang):
    psi4.energy('ccsd(t)')
    return [psi4.variable("SCF TOTAL ENERGY"), psi4.variable("MP2 TOTAL ENERGY"),psi4.variable("CCSD TOTAL ENERGY"),psi4.variable("CCSD(T) TOTAL ENERGY")]

def main(molname, memory, basis, bl_upper, bl_lower, bl_step, ang, num_threads):
    Path(molname).mkdir()
    psi4.set_output_file(f'{molname}/{molname}.psi4out', append=False)
    psi4.set_memory(f'{memory} MB')
    psi4.core.set_num_threads(num_threads)
    psi4.set_options({'basis':basis})
    num_points = int((bl_upper-bl_lower)/bl_step + 1)
    binding_data = np.zeros((num_points,6))
    i = 0
    for bl in np.linspace(bl_lower,bl_upper,num_points):
        dirname = f'{molname}/{bl:.2f}_{ang:.2f}'
        mol, wfn = generate_water(bl, ang)
        generate_dat_psi(dirname, mol, wfn)
        e_list = run_psi4(bl, ang)
        binding_data[i,:2] = [bl, ang]
        binding_data[i,2:] = e_list
        i += 1
        with open(f'{dirname}/reference.dat','w') as f:
            f.write(f'HF: {e_list[0]}\n')
            f.write(f'MP2: {e_list[1]}\n')
            f.write(f'CCSD: {e_list[2]}\n')
            f.write(f'CCSD(T): {e_list[3]}\n')
    np.savetxt(f'{molname}/binding_data.dat',binding_data,['%5.3f','%6.3f','%17.15f','%17.15f','%17.15f','%17.15f'])

if __name__ == '__main__':
    # Name of molecule
    molname = 'n2'
    
    # MBs
    memory = 2000

    # Threading
    num_threads = 8
    
    # Basis set to be used
    basis = 'cc-pvtz'
    molname = f'{molname}-{basis}'
    
    # Angstrom
    bl_upper = 2.0
    bl_lower = 2.0
    bl_step = 0.02
    
    # Degrees
    ang = 104.45
    
    main(molname, memory, basis, bl_upper, bl_lower, bl_step, ang, num_threads)
