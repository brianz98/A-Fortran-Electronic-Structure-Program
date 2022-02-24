#!/usr/bin/env python
# coding: utf-8

import numpy as np
from pathlib import Path
import subprocess as sp

def generate_gamess_input_file(bl, ang, dirname, calcname):
    geom_string = f"""
    cnv 2

    O
    H 1 {bl}
    H 1 {bl} 2 {ang}\n
    """

    with open(f'{dirname}/{calcname}.inp','w') as f:
        f.write(' $contrl scftyp=rhf coord=zmt runtyp=energy units=angs cctyp=cr-cc ispher=1 $end\n')
        f.write(' $system mwords=5 memddi=10 $end\n')
        f.write(' $guess  guess=huckel $end\n')
        f.write(' $basis  gbasis=ccd $end\n')
        f.write(' $data\n')
        f.write(geom_string)
        f.write(' $end')

def run_gamess(bl, ang, gamess_dir, directory, calc_name):
    
    generate_gamess_input_file(bl, ang, directory, calc_name)
    gamess_string = gamess_dir+f' {calc_name} 00 1 1 1'
    result = sp.check_output(gamess_string,cwd=directory, shell=True).decode('utf-8').split('\n')

    energy = np.zeros(11)
    with open(f'{directory}/{calc_name}.out', 'w') as f:
        for idx, line in enumerate(result):
            f.write(line+'\n')
            if 'REFERENCE ENERGY:' in line:
                energy[0] = float(line.split(' ')[-1])
            if 'MBPT(2) ENERGY:' in line:
                energy[1] = float(line.split('   CORR.E')[0].split(' ')[-1])
            if 'CCSD    ENERGY:' in line:
                energy[1] = float(line.split('   CORR.E')[0].split(' ')[-1])
            if 'CCSD[T] ENERGY:' in line:
                energy[2] = float(line.split('   CORR.E')[0].split(' ')[-1])
            if 'CCSD(T) ENERGY:' in line:
                energy[3] = float(line.split('   CORR.E')[0].split(' ')[-1])
            if 'R-CCSD[T] ENERGY:' in line:
                energy[4] = float(line.split('   CORR.E')[0].split(' ')[-1])
            if 'R-CCSD(T) ENERGY:' in line:
                energy[5] = float(line.split('   CORR.E')[0].split(' ')[-1])
            if 'CR-CCSD[T] ENERGY:' in line:
                energy[6] = float(line.split('   CORR.E')[0].split(' ')[-1])
            if 'CR-CCSD(T) ENERGY:' in line:
                energy[7] = float(line.split('   CORR.E')[0].split(' ')[-1])
            if 'T1 DIAGNOSTIC' in line:
                energy[8] = float(line.split(' ')[-1])
            if 'R-CCSD[T] DENOMINATOR' in line:
                energy[9] = float(line.split(' ')[-1])
            if 'R-CCSD(T) DENOMINATOR' in line:
                energy[10] = float(line.split(' ')[-1])
    return energy

def main(molname, basis, bl_upper, bl_lower, bl_step, ang, gamess_dir, gamess_scrdir):
    # Clear scratch
    sp.call(f'rm {gamess_scrdir}/*',shell=True)

    Path(molname).mkdir(exist_ok=True)

    num_points = int((bl_upper-bl_lower)/bl_step + 1)
    # [bondlength, angle, e_hf, e_mp2, e_ccsd, e_ccsd[t], e_ccsd(t), e_rccsd[t], e_rccsd(t), e_crccsd[t], e_crccsd(t), t1_diagn, [t]_denom, (t)_denom]
    binding_data = np.zeros((num_points,13))
    i = 0
    for bl in np.linspace(bl_lower,bl_upper,num_points):
        dirname = f'{molname}'
        calc_name = f'{molname}_{bl}_{ang}'
        e_list = run_gamess(bl, ang, gamess_dir, dirname, calc_name)
        binding_data[i,:2] = [bl, ang]
        binding_data[i,2:] = e_list
        i += 1

    fmt_arr = ['%5.3f','%6.3f'] + ['%17.15f']*11
    np.savetxt(f'{molname}/binding_data.dat',binding_data,fmt_arr)

if __name__ == '__main__':
    # Directory of the rungms script
    gamess_dir = '/scratch/zz376/software/gamess-stable/rungms'
    
    # Directory of GAMESS scratch, we need to clear it before each run
    gamess_scrdir = '/scratch/zz376/gamess/stable-userscr'

    # Name of molecule
    molname = 'h2'
        
    # Basis set to be used
    basis = 'sto-3g'
    molname = f'{molname}-{basis}'
    
    # Angstrom
    bl_upper = 1.0
    bl_lower = 0.6
    bl_step = 0.1
    
    # Degrees
    ang = 104.45
    
    main(molname, basis, bl_upper, bl_lower, bl_step, ang, gamess_dir, gamess_scrdir)
