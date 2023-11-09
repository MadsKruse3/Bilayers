from pathlib import Path
from utils import listdirs
import os
from asr.core import read_json
import matplotlib.pyplot as plt
import numpy as np
from ase.io import read
from gpaw import GPAW

def get_energy_difference(folder):
    energies = []
    name_energy = []
    atoms = read(f"{folder}/structure.json")
    n = len(atoms.get_chemical_symbols())
    dpath = f"{folder}/results-asr.bilayer_magnetism.json"
    if os.path.exists(dpath):
        data = read_json(dpath)
        
        AFM_norm = data["AFMNormIsZero"]
        print(AFM_norm)

        energydiff = data["eDIFF"] 
        if 'energydiff' in locals():
            energydiff = energydiff/n

        name_energy.append(("/".join((str(folder).split("/")[-2:])), abs(energydiff) ))
        energies.append(energydiff)
    
    return energies, name_energy, AFM_norm

def magnetic_type(folder):
    magmom_path = f"{folder}/gs_U3.0.gpw"
    if os.path.exists(magmom_path):
        calc = GPAW(f'{folder}/gs_U3.0.gpw', txt=None)
        magmom = calc.get_magnetic_moment()
        magmoms = calc.get_magnetic_moments()
        return magmoms, magmom
    else:
        magmom = []
        magmoms = []
        return magmoms, magmom

def get_exfoliation_energy(folder):
    energies = []
    name_energy = []
    atoms = read(f"{folder}/structure.json")
    dpath = f"{folder}/results-asr.bilayersummary.json"
    if os.path.exists(f"{folder}/results-asr.bilayersummary.json"):
        data = read_json(dpath)
        ex_energy = data["exfoliation_energy"] 

        name_energy.append(("/".join((str(folder).split("/")[-2:])), abs(ex_energy) ))
        energies.append(ex_energy)
    
    return energies, name_energy

def get_gap_type(folder):
    name_material = []
    atoms = read(f"{folder}/structure.json")
    dpath = f"{folder}/results-asr.gs.json"
    if os.path.exists(dpath):
        data = read_json(dpath)
        gap = data["gap"]
        if 'gap' in locals():
            name_material.append(("/".join((str(folder).split("/")[-2:])), gap ))
        else:
            gap = 'missing'
            name_material.append(("/".join((str(folder).split("/")[-2:])), gap ))
   
        return gap, name_material
    else:
        name_material.append(("/".join((str(folder).split("/")[-2:])), 'missing' ))
        return 'missing', name_material

def check_magmoms(folder):
    name_material = []
    atoms = read(f"{folder}/structure.json")
    dpath = f"{folder}/gs_U3.0.gpw"
    magnetic_atoms = []

    if os.path.exists(dpath):
    
        TM3d_atoms = ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']
        
        for atom in atoms:
            if atom.symbol in TM3d_atoms:
                magnetic_atoms.append(1)
            else:
                magnetic_atoms.append(0)
                
        magnetic_atoms_upper_layer = magnetic_atoms[:len(magnetic_atoms)//2]
        mag_atoms = [x for x, z in enumerate(magnetic_atoms) if z == 1]

        lines = []
        magmoms = []
        with open(f'{folder}/fm_U3.0.txt', 'r') as read_obj:
            for line in read_obj:
                line = line.strip()
                lines.append(line)
                if line == 'Local magnetic moments:': 
                    mag_line_fm = line
                    
            if not 'mag_line_fm' in locals():
                read_obj.close()
                deviation_matrix_fm = 'None'
                deviation_matrix_afm = 'None'
                return deviation_matrix_fm, deviation_matrix_afm

            start_line = lines.index(mag_line_fm)
            maglines = [start_line + 1 + x for x in mag_atoms]

            for line in lines:
                if lines.index(line) in maglines:
                    newline = str(str(line).split("(")[-1]).split(",")[-1]
                    magmoms.append(str(newline).split(")")[0])
    
            magmoms = [float(x) for x in magmoms]
            deviation_matrix_fm = []
            for x in magmoms:
                deviation_matrix_fm.append([ (abs(x) - abs(y))/abs(x) for y in magmoms])
        
            deviation_matrix_fm = np.array(deviation_matrix_fm)
            
        read_obj.close()

        lines = []
        magmoms = []
        with open(f'{folder}/afm_U3.0.txt', 'r') as read_obj:
            for line in read_obj:
                line = line.strip()
                lines.append(line)
                if line == 'Local magnetic moments:': 
                    mag_line_afm = line

            if not 'mag_line_afm' in locals():
                read_obj.close()
                deviation_matrix_fm = 'None'
                deviation_matrix_afm = 'None'
                return deviation_matrix_fm, deviation_matrix_afm
                
            start_line = lines.index(mag_line_afm) ##AFM
            maglines = [start_line + 1 + x for x in mag_atoms]

            for line in lines:
                if lines.index(line) in maglines:
                    newline = str(str(line).split("(")[-1]).split(",")[-1]
                    magmoms.append(str(newline).split(")")[0])

            magmoms = [float(x) for x in magmoms]
            deviation_matrix_afm = []
            for x in magmoms:
                deviation_matrix_afm.append([ (abs(x) - abs(y))/abs(x) for y in magmoms])
        
            deviation_matrix_afm = np.array(deviation_matrix_afm)
    
        read_obj.close()
    
        return deviation_matrix_fm, deviation_matrix_afm

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("folders", nargs="*", help="Monolayer folders to analyse.")

    args = parser.parse_args()
    
    if len(args.folders) > 0:
        folders = [Path(x).absolute() for x in args.folders]
    else:
        folders = [Path(".").absolute()]

    Ediffs = []
    name_energies = []
    mag_moms = []
    mag_mom = []
    ex_energy = []

    insulators_exfoliation = []
    insulators_exchange = []
    metals_exfoliation = []
    metals_exchange = []
    
    
    for folder in folders:
        if os.path.isfile(f'{folder}/fm_U3.0.txt') and os.path.isfile(f'{folder}/afm_U3.0.txt') and os.path.isfile(f'{folder}/gs_U3.0.gpw'):
            print(folder)
            Ediff, name_energy, AFM_norm = get_energy_difference(folder)
            magmoms, magmom = magnetic_type(folder)
            gap, name_material = get_gap_type(folder)
            deviation_matrix_fm, deviation_matrix_afm = check_magmoms(folder)

            #if abs(Ediff[0]) > 0.025:
            #    print(name_energy)

            check_values = []
            if not deviation_matrix_fm == 'None':
                if not deviation_matrix_afm == 'None':
                    for x in deviation_matrix_fm, deviation_matrix_afm:
                        for y in x:
                            for z in y:
                                if abs(z) > 0.05:
                                    check_values.append(z)

            if deviation_matrix_fm == 'None' or deviation_matrix_afm == 'None':
                check_values.append('None')

            if len(check_values) == 0:
                if AFM_norm == True:
                    exfoliation_energy, exfoliation_energy_material = get_exfoliation_energy(folder)
                    ex_energy.extend(exfoliation_energy)

                    if gap == 0.0:
                        metals_exchange.append(Ediff)
                        metals_exfoliation.append(exfoliation_energy)
                    else:
                        insulators_exchange.append(Ediff)
                        insulators_exfoliation.append(exfoliation_energy)

                    Ediffs.extend(Ediff)
                    name_energies.extend(name_energy)
        

    os.chdir('/home/niflheim2/cmr/WIP/stacking/tree-mads/data_and_plots/')
    np.savetxt('CrI3_AA_insulators_exfoliation_energy.npy', insulators_exfoliation)   
    np.savetxt('CrI3_AA_insulators_exchange_energy.npy', insulators_exchange)
    np.savetxt('CrI3_AA_metals_exfoliation_energy.npy', metals_exfoliation)   
    np.savetxt('CrI3_AA_metals_exchange_energy.npy', metals_exchange)

    
