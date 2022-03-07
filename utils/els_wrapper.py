#!/usr/bin/env python
# coding: utf-8

import psi4
import numpy as np
from pathlib import Path
import subprocess as sp

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

    sp.call(f'cp els.in {dirname}',shell=True)

def generate_water(bl, ang):
    mol = psi4.geometry(f""" 
        F
        F 1 {bl} 
        units = angstrom
        """)
    wfn = psi4.core.Wavefunction(mol, psi4.core.BasisSet.build(mol))
    return mol, wfn

def run_psi4(bl, ang):
    psi4.energy('ccsd(t)')
    return [psi4.variable("SCF TOTAL ENERGY"), psi4.variable("MP2 TOTAL ENERGY"),psi4.variable("CCSD TOTAL ENERGY"),psi4.variable("CCSD(T) TOTAL ENERGY")]

def run_els(els_dir, directory):
    result = sp.check_output(els_dir,cwd=directory).decode('utf-8').split('\n')
    energy = np.zeros(12)
    with open(f'{directory}/els.out', 'w') as f:
        for idx, line in enumerate(result):
            f.write(line+'\n')
            if 'RHF energy:' in line:
                energy[0] = float(line.split(' ')[-1])
            if 'MP2 energy:' in line:
                energy[1] = float(line.split(' ')[-1])
            if ' CCSD energy:' in line:
                energy[2] = float(line.split(' ')[-1])
            if ' CCSD[T] energy:' in line:
                energy[3] = float(line.split(' ')[-1])
            if ' CCSD(T) energy:' in line:
                energy[4] = float(line.split(' ')[-1])
            if ' R-CCSD[T] energy:' in line:
                energy[5] = float(line.split(' ')[-1])
            if ' R-CCSD(T) energy:' in line:
                energy[6] = float(line.split(' ')[-1])
            if ' CR-CCSD[T] energy:' in line:
                energy[7] = float(line.split(' ')[-1])
            if ' CR-CCSD(T) energy:' in line:
                energy[8] = float(line.split(' ')[-1])
            if ' T1 diagnostic:' in line:
                energy[9] = float(line.split(' ')[-1])
            if ' D[T]:' in line:
                energy[10] = float(line.split(' ')[-1])
            if ' D(T):' in line:
                energy[11] = float(line.split(' ')[-1])
    return energy

def main(molname, memory, basis, bl_upper, bl_lower, bl_step, ang, els_dir):
    Path(molname).mkdir(exist_ok=True)
    psi4.set_output_file(f'{molname}/{molname}.psi4out', append=False)
    psi4.set_memory(f'{memory} MB')
    psi4.set_options({'basis':basis})
    num_points = round((bl_upper-bl_lower)/bl_step + 1)
    binding_data_psi4 = np.zeros((num_points,6))
    binding_data_els = np.zeros((num_points,14))
    i = 0
    for bl in np.linspace(bl_lower,bl_upper,num_points):
        dirname = f'{molname}/{bl:.2f}_{ang:.2f}'
        mol, wfn = generate_water(bl, ang)
        generate_dat_psi(dirname, mol, wfn)

        e_list = run_psi4(bl, ang)
        binding_data_psi4[i,:2] = [bl, ang]
        binding_data_psi4[i,2:] = e_list
        with open(f'{dirname}/reference.dat','w') as f:
            f.write(f'HF: {e_list[0]}\n')
            f.write(f'MP2: {e_list[1]}\n')
            f.write(f'CCSD: {e_list[2]}\n')
            f.write(f'CCSD(T): {e_list[3]}\n')

        e_list_els = run_els(els_dir, dirname)
        binding_data_els[i,:2] = [bl, ang]
        binding_data_els[i,2:] = e_list_els
        with open(f'{dirname}/els_energy.dat','w') as f:
            f.write(f'HF: {e_list_els[0]}\n')
            f.write(f'MP2: {e_list_els[1]}\n')
            f.write(f'CCSD: {e_list_els[2]}\n')
            f.write(f'CCSD[T]: {e_list_els[3]}\n')
            f.write(f'CCSD(T): {e_list_els[4]}\n')
            f.write(f'R-CCSD[T]: {e_list_els[5]}\n')
            f.write(f'R-CCSD(T): {e_list_els[6]}\n')
            f.write(f'CR-CCSD[T]: {e_list_els[7]}\n')
            f.write(f'CR-CCSD(T): {e_list_els[8]}\n')
            f.write(f'T1 diagnostic: {e_list_els[9]}\n')
            f.write(f'D[T]: {e_list_els[10]}\n')
            f.write(f'D(T): {e_list_els[11]}\n')
        i += 1

    np.savetxt(f'{molname}/binding_data_psi4.dat',binding_data_psi4,['%5.3f','%6.3f','%17.15f','%17.15f','%17.15f','%17.15f'])
    fmt_arr = ['%5.3f','%6.3f'] + ['%17.15f']*12
    np.savetxt(f'{molname}/binding_data_els.dat',binding_data_els,fmt_arr)

if __name__ == '__main__':
    # Directory of els.x
    els_dir = '/mnt/c/Users/zdj51/Documents/code/electronic-structure/els.x'

    # Name of molecule
    molname = 'f2'
    
    # MBs
    memory = 2000
    
    # Basis set to be used
    basis = 'cc-pvdz'
    molname = f'{molname}-{basis}'
    
    # Angstrom
    bl_lower = 1.0
    bl_upper = 3.0
    bl_step = 0.05
    
    # Degrees
    ang = 0
    
    main(molname, memory, basis, bl_upper, bl_lower, bl_step, ang, els_dir)
