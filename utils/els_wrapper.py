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

def generate_water(bl, ang):
    mol = psi4.geometry(f""" 
        O
        H 1 {bl}
        H 1 {bl} 2 {ang}
        units = angstrom
        """)
    wfn = psi4.core.Wavefunction(mol, psi4.core.BasisSet.build(mol))
    return mol, wfn

def run_psi4(bl, ang, read_in):
    if not read_in:
        e_ccsd_t, wfn = psi4.energy('ccsd(t)', return_wfn=True)
        wfn.to_file('psi4_out')
    else:
        sp.call(f'mv psi4_out.npy psi4_in.npy', shell=True)
        e_ccsd_t, wfn = psi4.energy('ccsd(t)',restart_file='psi4_in',return_wfn=True)
        wfn.to_file('psi4_out')

    psi4.core.clean()
    return [psi4.variable("SCF TOTAL ENERGY"), psi4.variable("MP2 TOTAL ENERGY"),psi4.variable("CCSD TOTAL ENERGY"),psi4.variable("CCSD(T) TOTAL ENERGY")]

def run_els(els_dir, directory, previous_directory, read_in):
    if not read_in:
        sp.call(f'cp els_noread.in {directory}/els.in',shell=True)
    else:
        sp.call(f'cp els.in {directory}',shell=True)
        sp.call(f'cp {previous_directory}/guess_out.dat {directory}/guess_in.dat',shell=True)

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

def main(molname, memory, basis, bl_upper, bl_lower, bl_step, ang, els_dir, read_in, num_threads, psi4_off, psi4_read_in):
    Path(molname).mkdir(exist_ok=True)
    psi4.set_output_file(f'{molname}/{molname}.psi4out', append=False)
    psi4.set_memory(f'{memory} MB')
    psi4.core.set_num_threads(num_threads)
    psi4.set_options({'basis':basis,'ccenergy__maxiter':200,'scf__maxiter':200})
    num_points = round((bl_upper-bl_lower)/bl_step + 1)
    binding_data_psi4 = np.zeros((num_points,6))
    binding_data_els = np.zeros((num_points,14))
    

    i = 0
    # We need this for read in/out of SCF guesses
    prev_dirname = ''

    for bl in np.linspace(bl_lower,bl_upper,num_points):
        dirname = f'{molname}/{bl:.2f}_{ang:.2f}'
        print(f'Doing calculations in {dirname}')
        
        mol, wfn = generate_water(bl, ang)
        generate_dat_psi(dirname, mol, wfn)

        try:
            if i == 0:
                e_list = run_psi4(bl, ang, False)
            else:
                e_list = run_psi4(bl, ang, psi4_read_in)

        except:
            if psi4_off:
                # Don't care, continue
                print('Psi4 failure, we\'ll stop running reference calculations')
                continue
            else:
                print('Psi4 failure, stopping')
                break

        binding_data_psi4[i,:2] = [bl, ang]
        binding_data_psi4[i,2:] = e_list
        with open(f'{dirname}/reference.dat','w') as f:
            f.write(f'HF: {e_list[0]}\n')
            f.write(f'MP2: {e_list[1]}\n')
            f.write(f'CCSD: {e_list[2]}\n')
            f.write(f'CCSD(T): {e_list[3]}\n')

        try:
            if i == 0: 
                # Must set scf_read_in setting to .false. since there's no previous calculation
                e_list_els = run_els(els_dir, dirname, prev_dirname, read_in=False)
            else:
                # If read_in not specified, never read in
                e_list_els = run_els(els_dir, dirname, prev_dirname, read_in=read_in)
        except:
            print('els failure')
            break

        prev_dirname = dirname
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
    molname = 'h2o'
    
    # MBs
    memory = 2000

    # Parallelisation in Psi4
    num_threads = 8
    
    # Basis set to be used
    basis = 'cc-pvdz'
    molname = f'{molname}-{basis}'
    
    # Angstrom
    bl_lower = 2.00
    bl_upper = 2.2
    bl_step = 0.02
    
    # Degrees
    ang = 104.45

    # Read in guesses?
    read_in = True

    # Psi4 being trash and not converging? Don't run reference calculations
    psi4_off = False

    # Also do read-in for Psi4?
    psi4_read_in = True
    
    main(molname, memory, basis, bl_upper, bl_lower, bl_step, ang, els_dir, read_in, num_threads, psi4_off, psi4_read_in)
