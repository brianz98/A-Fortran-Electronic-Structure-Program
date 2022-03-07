#!/usr/bin/env python
# coding: utf-8

import numpy as np
from pathlib import Path
import subprocess as sp

def generate_gamess_input_file(bl, ang, dirname, calc_name, basis):
    geom_string = f"""
dnh 2

F
F 1 {bl}\n"""

    with open(f'{dirname}/{calc_name}.inp','w') as f:
        f.write(' $contrl scftyp=rhf coord=zmt runtyp=energy units=angs cctyp=cr-cc ispher=1 $end\n')
        f.write(' $system mwords=100 memddi=500 $end\n')
        f.write(' $guess  guess=huckel $end\n')
        f.write(' $ccinp  maxcc=100 ncore=0 $end\n')
        f.write(f' $basis  {basis} $end\n')
        f.write(' $data\n')
        f.write(geom_string)
        f.write(' $end')

def run_gamess(bl, ang, gamess_dir, directory, calc_name, basis):
    
    generate_gamess_input_file(bl, ang, directory, calc_name, basis)
    gamess_string = gamess_dir+f' {calc_name} 00 1 1 1'
    result = sp.check_output(gamess_string,cwd=directory, shell=True).decode('utf-8').split('\n')

    energy = np.zeros(12)
    with open(f'{directory}/{calc_name}.out', 'w') as f:
        for idx, line in enumerate(result):
            f.write(line+'\n')
            # Some spaces applied for the case the RCCSD(T) gets recognised as CCSD(T), for example.
            if 'REFERENCE ENERGY:' in line:
                energy[0] = float(line.split(' ')[-1])
            if 'MBPT(2) ENERGY:' in line:
                energy[1] = float(line.split('   CORR.E')[0].split(' ')[-1])
            if ' CCSD    ENERGY:' in line:
                energy[2] = float(line.split('   CORR.E')[0].split(' ')[-1])
            if ' CCSD[T] ENERGY:' in line:
                energy[3] = float(line.split('   CORR.E')[0].split(' ')[-1])
            if ' CCSD(T) ENERGY:' in line:
                energy[4] = float(line.split('   CORR.E')[0].split(' ')[-1])
            if ' R-CCSD[T] ENERGY:' in line:
                energy[5] = float(line.split('   CORR.E')[0].split(' ')[-1])
            if ' R-CCSD(T) ENERGY:' in line:
                energy[6] = float(line.split('   CORR.E')[0].split(' ')[-1])
            if 'CR-CCSD[T] ENERGY:' in line:
                energy[7] = float(line.split('   CORR.E')[0].split(' ')[-1])
            if 'CR-CCSD(T) ENERGY:' in line:
                energy[8] = float(line.split('   CORR.E')[0].split(' ')[-1])
            if 'T1 DIAGNOSTIC' in line:
                energy[9] = float(line.split(' ')[-1])
            if ' R-CCSD[T] DENOMINATOR' in line:
                energy[10] = float(line.split(' ')[-1])
            if ' R-CCSD(T) DENOMINATOR' in line:
                energy[11] = float(line.split(' ')[-1])
    return energy

def main(molname, basis, bl_upper, bl_lower, bl_step, ang_upper, ang_lower, ang_step, gamess_dir, gamess_scrdir):
    # Clear scratch
    sp.call(f'rm {gamess_scrdir}/*',shell=True)

    Path(molname).mkdir(exist_ok=True)

    num_bl_points = round((bl_upper-bl_lower)/bl_step + 1)
    num_ang_points = round((ang_upper-ang_lower)/ang_step + 1)
    # [bondlength, angle, e_hf, e_mp2, e_ccsd, e_ccsd[t], e_ccsd(t), e_rccsd[t], e_rccsd(t), e_crccsd[t], e_crccsd(t), t1_diagn, [t]_denom, (t)_denom]
    binding_data = np.zeros((num_bl_points*num_ang_points,14))
    i = 0
    if (ang_upper == ang_lower):
        ang = 0
        # Disable looping over angles
        for bl in np.linspace(bl_lower,bl_upper,num_bl_points):
            dirname = f'{molname}'
            calc_name = f'{molname}_{bl:.3f}'
            e_list = run_gamess(bl, ang, gamess_dir, dirname, calc_name, basis)
            binding_data[i,:2] = [bl, ang]
            binding_data[i,2:] = e_list
            i += 1
    else:
        # Loop over angles
        for bl in np.linspace(bl_lower,bl_upper,num_bl_points):
            for ang in np.linspace(ang_lower,ang_upper,num_ang_points):
                dirname = f'{molname}'
                calc_name = f'{molname}_{bl:.3f}_{ang:.3f}'
                e_list = run_gamess(bl, ang, gamess_dir, dirname, calc_name, basis)
                binding_data[i,:2] = [bl, ang]
                binding_data[i,2:] = e_list
                i += 1

    fmt_arr = ['%5.3f','%6.3f'] + ['%17.15f']*12
    np.savetxt(f'{molname}/binding_data.dat',binding_data,fmt_arr)

if __name__ == '__main__':
    # Directory of the rungms script
    gamess_dir = '/scratch/zz376/software/gamess-stable/rungms'
    
    # Directory of GAMESS scratch, we need to clear it before each run
    gamess_scrdir = '/scratch/zz376/gamess/stable-userscr'

    # Name of molecule
    molname = 'f2'
        
    # Angstrom
    bl_lower = 1.0
    bl_upper = 3.0
    bl_step = 0.02
    
    # Degrees
    ang_upper = 150
    ang_lower = 150
    ang_step = 2.5
    
    # Basis set to be used
    basis_dict = {
        'sto-3g':'gbasis=sto ngauss=3', 
        '3-21g':'gbasis=n21 ngauss=3',
        'cc-pvdz':'gbasis=ccd',
        'cc-pvtz':'gbasis=cct',
        'cc-pvqz':'gbasis=ccq',
        'aug-cc-pvdz':'gbasis=accd',
        'aug-cc-pvtz':'gbasis=acct',
        'aug-cc-pvqz':'gbasis=accq',
        'cc-pcvdz':'gbasis=ccdc',
        'cc-pcvtz':'gbasis=cctc',
        'cc-pcvqz':'gbasis=ccqc',
        'aug-cc-pcvdz':'gbasis=accdc',
        'aug-cc-pcvtz':'gbasis=acctc',
        'aug-cc-pcvqz':'gbasis=accqc',
        'cc-pwcvdz':'gbasis=ccdwc',
        'cc-pwcvtz':'gbasis=cctwc',
        'cc-pwcvqz':'gbasis=ccqwc',
        'aug-cc-pwcvdz':'gbasis=accdwc',
        'aug-cc-pwcvtz':'gbasis=acctwc',
        'aug-cc-pwcvqz':'gbasis=accqwc'
    }
    
    for basis in ['cc-pvdz','aug-cc-pvqz']:
        dirname = f'{molname}-{basis}'
        gamess_basis_string = basis_dict[basis]
        
        main(dirname, gamess_basis_string, bl_upper, bl_lower, bl_step, ang_upper, ang_lower, ang_step, gamess_dir, gamess_scrdir)
