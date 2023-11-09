from pathlib import Path
#from utils import listdirs
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
    dpath = f"{folder}/results-asr.interlayer_magnetic_exchange.json"
    if os.path.exists(dpath):
        data = read_json(dpath)

        energydiff = data["eDIFF"] 
        if 'energydiff' in locals():
            energydiff = energydiff/n

        name_energy.append(("/".join((str(folder).split("/")[-2:])), abs(energydiff) ))
        energies.append(energydiff)
        
    return energies, name_energy

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
    dpath = f"{folder}/results-asr.interlayer_magnetic_exchange.json"
    magnetic_atoms = []

    if os.path.exists(dpath):
    
        TM3d_atoms = ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']
        
        for atom in atoms:
            if atom.symbol in TM3d_atoms:
                magnetic_atoms.append(1)
            else:
                magnetic_atoms.append(0)
                
        mag_atoms = [x for x, z in enumerate(magnetic_atoms) if z == 1]

        data = read_json(dpath)

        magmoms_FM = data["magmoms_fm"]
        magmom_FM = data["M_FM"]
        magmoms_AFM = data["magmoms_afm"]
        magmom_AFM = data["M_AFM"]
 
        lines = []
        magmoms = []
        
        deviation_matrix_fm = []
        for x in mag_atoms:
            deviation_matrix_fm.append([ (abs(magmoms_FM[x]) - abs(magmoms_FM[y]))/abs(magmoms_FM[x]) for y in mag_atoms])
        deviation_matrix_fm = np.array(deviation_matrix_fm)
        
        deviation_matrix_afm = []
        for x in mag_atoms:
            deviation_matrix_afm.append([ (abs(magmoms_AFM[x]) - abs(magmoms_AFM[y]))/abs(magmoms_AFM[x]) for y in mag_atoms])
        deviation_matrix_afm = np.array(deviation_matrix_afm)
            
        if not np.allclose((magmom_AFM/len(mag_atoms)),0, atol=0.01):
            deviation_matrix_fm = 'missing'
            deviation_matrix_afm = 'missing'
   
        #FM_syms = bilayer_FM.get_chemical_symbols()
        FM_sign = []
        negative_count = []
        for x in mag_atoms:
            FM_sign.append(np.sign(magmoms_FM[x]))
                
        for x in FM_sign:
            if x == -1:
                negative_count.append(1)

        if not len(negative_count) == 0:
            deviation_matrix_fm = 'missing'
            deviation_matrix_afm = 'missing'

    else:
        deviation_matrix_fm = 'missing'
        deviation_matrix_afm = 'missing'


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
    
    outliers = []
    FM = []
    AFM = []
    
    FM_metals = []
    FM_insulators = []
    AFM_metals = []
    AFM_insulators = []

    for folder in folders:
        if os.path.isfile(f'{folder}/results-asr.interlayer_magnetic_exchange.json'):
            Ediff, name_energy = get_energy_difference(folder)
            gap, name_material = get_gap_type(folder)
            deviation_matrix_fm, deviation_matrix_afm = check_magmoms(folder)

            check_values = []
            if not deviation_matrix_fm == 'missing':
                if not deviation_matrix_afm == 'missing':
                    for x in deviation_matrix_fm, deviation_matrix_afm:
                        for y in x:
                            for z in y:
                                if abs(z) > 0.05:
                                    check_values.append(z)

            if deviation_matrix_fm == 'missing' or deviation_matrix_afm == 'missing':
                check_values.append('missing')


            if len(check_values) == 0:

                exfoliation_energy, exfoliation_energy_material = get_exfoliation_energy(folder)
                ex_energy.extend(exfoliation_energy)
                
                if gap == 0.0:
                    metals_exchange.append(Ediff)
                    metals_exfoliation.append(exfoliation_energy)
                    if Ediff[0] < 0:
                        FM_metals.append(Ediff)
                    if Ediff[0] > 0:
                        AFM_metals.append(Ediff)
                else:
                    insulators_exchange.append(Ediff)
                    insulators_exfoliation.append(exfoliation_energy)
                    if Ediff[0] < 0:
                        FM_insulators.append(Ediff[0])
                    if Ediff[0] > 0:
                        AFM_insulators.append(Ediff[0])
                    
                Ediffs.extend(Ediff)
                name_energies.extend(name_energy)
                
                #print(Ediff[0]*1e3)
                #if 1e3*abs(Ediff[0]) < 5:
                    #outliers.append(name_energy)
                #    print(name_energy)

                if Ediff[0] > 0:
                    AFM.append(Ediff[0])
                if Ediff[0] < 0: 
                    FM.append(Ediff[0])

    
    os.chdir('/home/niflheim2/cmr/WIP/stacking/tree-mads/data_and_plots/')
    np.savetxt('insulators_exfoliation_energy.npy', insulators_exfoliation)   
    np.savetxt('insulators_exchange_energy.npy', insulators_exchange)
    np.savetxt('metals_exfoliation_energy.npy', metals_exfoliation)   
    np.savetxt('metals_exchange_energy.npy', metals_exchange)
    #np.savetxt('outliers_exchange_energy.npy', outliers)

    print('Number of metals:', len(metals_exchange))
    print('Number of insulators:', len(insulators_exchange))
    #print('Outlier materials:', outliers)
    print('FM bilayers:', len(FM))
    print('AFM bilayers:', len(AFM))

    print('FM insulators', len(FM_insulators))
    print('FM metals', len(FM_metals))
    print('AFM insulators', len(AFM_insulators))
    print('AFM metals', len(AFM_metals))
