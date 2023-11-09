from pathlib import Path
from utils import listdirs
import os, re
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

def get_ehull(monolayer_folder):
    fname = 'results-asr.convex_hull.json'
    path = "/".join([x for x in monolayer_folder.split("/") if x != ""][-3:])
    c2db_path = "/home/niflheim2/cmr/C2DB-ASR/tree/" + path + "/" + fname
    icsd_cod_path = "/home/niflheim2/cmr/C2DB-ASR/icsd_cod_materials/tree/" + path + "/" + fname
    excluded_AFM_path = "/home/niflheim2/cmr/C2DB-ASR/excluded_AFM_tree/" + path + "/" + fname
    if Path(c2db_path).is_file():
        data = read_json(f'{c2db_path}')
        ehull = data["ehull"]
    
    if Path(icsd_cod_path).is_file():
        data = read_json(f'{icsd_cod_path}')
        ehull = data["ehull"]
    
    if Path(excluded_AFM_path).is_file():
        data = read_json(f'{excluded_AFM_path}')
        ehull = data["ehull"]

    if not Path(c2db_path).is_file() and not Path(icsd_cod_path).is_file() and not Path(excluded_AFM_path).is_file():
        print(monolayer_folder)
        ehull = 'missing'
    
    return ehull

def get_IDs(monolayer_folder):
    fname = 'info.json'
    path = "/".join([x for x in monolayer_folder.split("/") if x != ""][-3:])
    c2db_path = "/home/niflheim2/cmr/C2DB-ASR/tree/" + path + "/" + fname
    icsd_cod_path = "/home/niflheim2/cmr/C2DB-ASR/icsd_cod_materials/tree/" + path + "/" + fname
    excluded_AFM_path = "/home/niflheim2/cmr/C2DB-ASR/excluded_AFM_tree/" + path + "/" + fname
    if Path(c2db_path).is_file():
        data = read_json(f'{c2db_path}')
        try:
            ID_icsd = data["icsd_id"]
        except Exception:
            pass
        try:
            ID_cod = data["cod_id"]
        except Exception:
            pass
        if not 'ID_cod' in locals():
            ID_cod = 'missing'
        if not 'ID_icsd' in locals():
            ID_icsd = 'missing'
    
    if Path(icsd_cod_path).is_file():
        data = read_json(f'{icsd_cod_path}')
        try:
            ID_icsd = data["icsd_id"]
        except Exception:
            pass
        try:
            ID_cod = data["cod_id"]
        except Exception:
            pass
        if not 'ID_cod' in locals():
            ID_cod = 'missing'
        if not 'ID_icsd' in locals():
            ID_icsd = 'missing'
    
    if Path(excluded_AFM_path).is_file():
        data = read_json(f'{excluded_AFM_path}')
        try:
            ID_icsd = data["icsd_id"]
        except Exception:
            pass
        try:
            ID_cod = data["cod_id"]
        except Exception:
            pass
        if not 'ID_cod' in locals():
            ID_cod = 'missing'
        if not 'ID_icsd' in locals():
            ID_icsd = 'missing'

    if not Path(c2db_path).is_file() and not Path(icsd_cod_path).is_file() and not Path(excluded_AFM_path).is_file():
        ID_icsd = 'missing'
        ID_cod = 'missing'
    
    return ID_icsd, ID_cod        
    
def get_monolayer_id(monolayer_folder):
    monolayer_id = read_json(f"{monolayer_folder}/structure.json")
    monolayer_id_data = monolayer_id[1]
    data_id = monolayer_id_data["unique_id"]
    return data_id

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
            
        if np.allclose((magmom_AFM/len(mag_atoms)),0, atol=0.01):
            deviation_matrix_fm == 'missing'
            deviation_matrix_afm == 'missing'
   
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
    ex_energy = []

    energies = []

    unique_monolayers = []
    monolayers = []
    
    all_FM_bilayers = []
    all_AFM_bilayers = []
    mixed_bilayers = []
    
    monolayer_bulk_origin_bilayers_all_FM = []
    monolayer_no_bulk_origin_but_small_ehull_bilayers_all_FM = []
    monolayer_remaining_bilayers_all_FM = []

    monolayer_bulk_origin_bilayers_all_AFM = []
    monolayer_no_bulk_origin_but_small_ehull_bilayers_all_AFM = []
    monolayer_remaining_bilayers_all_AFM = []

    monolayer_bulk_origin_bilayers_mixed = []
    monolayer_no_bulk_origin_but_small_ehull_bilayers_mixed = []
    monolayer_remaining_bilayers_mixed = []

    for folder in folders:
        if os.path.isfile(f'{folder}/results-asr.interlayer_magnetic_exchange.json'):

            monolayer_folder = re.split('[/]',f'{folder}')
            monolayer_folder.remove(monolayer_folder[-1])
            monolayer_folder = "/".join(monolayer_folder)

            Ediff, name_energy = get_energy_difference(folder)
            deviation_matrix_fm, deviation_matrix_afm = check_magmoms(folder)
            ID = get_monolayer_id(monolayer_folder)
                
            ID_icsd, ID_cod = get_IDs(monolayer_folder)

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

                #exfoliation_energy, exfoliation_energy_material = get_exfoliation_energy(folder)
                #ex_energy.extend(exfoliation_energy)

                if ID not in unique_monolayers:
                    unique_monolayers.append(ID)
                    if not len(energies) == 0:
                
                        unique_monolayers.append(ID)
                        if all([x < 0  for x in energies]) == True:
                            all_FM_bilayers.append(ID)
                            if not ID_icsd == 'missing' or not ID_cod == 'missing':
                                if not ID_icsd == 'missing':
                                    monolayer_bulk_origin_bilayers_all_FM.append(ID_icsd)
                                else:
                                    monolayer_bulk_origin_bilayers_all_FM.append(ID_cod)
                            else:
                                ehull = get_ehull(monolayer_folder)
                                if  type(ehull) == float: 
                                    if ehull < 50:
                                        monolayer_no_bulk_origin_but_small_ehull_bilayers_all_FM.append('missing')
                                else:
                                    monolayer_remaining_bilayers_all_FM.append('missing')

                        if all([x > 0  for x in energies]) == True:
                            all_AFM_bilayers.append(ID)
                            if not ID_icsd == 'missing' or not ID_cod == 'missing':
                                if not ID_icsd == 'missing':
                                    monolayer_bulk_origin_bilayers_all_AFM.append(ID_icsd)
                                else:
                                    monolayer_bulk_origin_bilayers_all_AFM.append(ID_cod)
                            else:
                                ehull = get_ehull(monolayer_folder)
                                if  type(ehull) == float: 
                                    if ehull < 50:
                                        monolayer_no_bulk_origin_but_small_ehull_bilayers_all_AFM.append('missing')
                                else:
                                    monolayer_remaining_bilayers_all_AFM.append('missing')

                        if not all([x < 0  for x in energies]) == True and not all([x > 0  for x in energies]) == True:
                            mixed_bilayers.append(ID)
                            if not ID_icsd == 'missing' or not ID_cod == 'missing':
                                if not ID_icsd == 'missing':
                                    monolayer_bulk_origin_bilayers_mixed.append(ID_icsd)
                                else:
                                    monolayer_bulk_origin_bilayers_mixed.append(ID_cod)
                            else:
                                ehull = get_ehull(monolayer_folder)
                                if  type(ehull) == float: 
                                    if ehull < 50:
                                        monolayer_no_bulk_origin_but_small_ehull_bilayers_mixed.append('missing')
                                else:
                                    monolayer_remaining_bilayers_mixed.append('missing')

                    energies = []
                    energies.append(Ediff[0])

                else:
                    energies.append(Ediff[0])

    if not len(energies) == 0:
        
        if all([x < 0  for x in energies]) == True:
            all_FM_bilayers.append(ID)
            if not ID_icsd == 'missing' or not ID_cod == 'missing':
                if not ID_icsd == 'missing':
                    monolayer_bulk_origin_bilayers_all_FM.append(ID_icsd)
                else:
                    monolayer_bulk_origin_bilayers_all_FM.append(ID_cod)
            else:
                ehull = get_ehull(monolayer_folder)
                if  type(ehull) == float: 
                    if ehull < 50:
                        monolayer_no_bulk_origin_but_small_ehull_bilayers_all_FM.append('missing')
                else:
                    monolayer_remaining_bilayers_all_FM.append('missing')

        if all([x > 0  for x in energies]) == True:
            all_AFM_bilayers.append(ID)
            if not ID_icsd == 'missing' or not ID_cod == 'missing':
                if not ID_icsd == 'missing':
                    monolayer_bulk_origin_bilayers_all_AFM.append(ID_icsd)
                else:
                    monolayer_bulk_origin_bilayers_all_AFM.append(ID_cod)
            else:
                ehull = get_ehull(monolayer_folder)
                if  type(ehull) == float: 
                    if ehull < 50:
                        monolayer_no_bulk_origin_but_small_ehull_bilayers_all_AFM.append('missing')
                else:
                    monolayer_remaining_bilayers_all_AFM.append('missing')

        if not all([x < 0  for x in energies]) == True and not all([x > 0  for x in energies]) == True:
            mixed_bilayers.append(ID)
            if not ID_icsd == 'missing' or not ID_cod == 'missing':
                if not ID_icsd == 'missing':
                    monolayer_bulk_origin_bilayers_mixed.append(ID_icsd)
                else:
                    monolayer_bulk_origin_bilayers_mixed.append(ID_cod)
            else:
                ehull = get_ehull(monolayer_folder)
                if  type(ehull) == float: 
                    if ehull < 50:
                        monolayer_no_bulk_origin_but_small_ehull_bilayers_mixed.append('missing')
                else:
                    monolayer_remaining_bilayers_mixed.append('missing')


    print('monolayers with purely FM bilayers:', len(all_FM_bilayers))
    print('monolayers with purely AFM bilayers:',len(all_AFM_bilayers))
    print('monolayers with mixed magnetic bilayers:',len(mixed_bilayers))
    
    print('monolayers with a bulk origin and purely FM bilayers:', len(monolayer_bulk_origin_bilayers_all_FM))
    print('monolayers without a bulk origin, but small ehull and purely FM bilayers:', len(monolayer_no_bulk_origin_but_small_ehull_bilayers_all_FM))
    print('monolayers remaining with purely FM bilayers:', len(monolayer_remaining_bilayers_all_FM))
    
    print('monolayers with a bulk origin and purely AFM bilayers:', len(monolayer_bulk_origin_bilayers_all_AFM))
    print('monolayers without a bulk origin, but small ehull and purely AFM bilayers:', len(monolayer_no_bulk_origin_but_small_ehull_bilayers_all_AFM))
    print('monolayers remaining with purely AFM bilayers:', len(monolayer_remaining_bilayers_all_AFM))

    print('monolayers with a bulk origin and mixed magnetic bilayers:', len(monolayer_bulk_origin_bilayers_mixed))
    print('monolayers without a bulk origin, but small ehull and mixed magnetic bilayers:', len(monolayer_no_bulk_origin_but_small_ehull_bilayers_mixed))
    print('monolayers remaining with mixed magnetic bilayers:', len(monolayer_remaining_bilayers_mixed))
  
    os.chdir('/home/niflheim2/cmr/WIP/stacking/tree-mads/data_and_plots/')

    np.savetxt('monolayers_with_all_FM_magnetic_bilayers.npy', all_FM_bilayers, fmt="%s")
    np.savetxt('monolayers_all_AFM_magnetic_bilayers.npy', all_AFM_bilayers, fmt="%s")
    np.savetxt('monolayers_mixed_magnetic_bilayers.npy', mixed_bilayers, fmt="%s")
    
    np.savetxt('monolayers_with_a_bulk_origin_and_all_FM_bilayers.npy', monolayer_bulk_origin_bilayers_all_FM, fmt="%s")
    np.savetxt('monolayers_without_a_bulk_origin_but_small_ehull_and_all_FM_bilayers.npy', monolayer_no_bulk_origin_but_small_ehull_bilayers_all_FM, fmt="%s")
    np.savetxt('monolayers_remaining_with_all_FM_bilayers.npy', monolayer_remaining_bilayers_all_FM, fmt="%s")

    np.savetxt('monolayers_with_a_bulk_origin_and_all_AFM_bilayers.npy', monolayer_bulk_origin_bilayers_all_AFM, fmt="%s")
    np.savetxt('monolayers_without_a_bulk_origin_but_small_ehull_and_all_AFM_bilayers.npy', monolayer_no_bulk_origin_but_small_ehull_bilayers_all_AFM, fmt="%s")
    np.savetxt('monolayers_remaining_with_all_AFM_bilayers.npy', monolayer_remaining_bilayers_all_AFM, fmt="%s")

    np.savetxt('monolayers_with_a_bulk_origin_and_mixed_bilayers.npy', monolayer_bulk_origin_bilayers_mixed, fmt="%s")
    np.savetxt('monolayers_without_a_bulk_origin_but_small_ehull_and_mixed_bilayers.npy', monolayer_no_bulk_origin_but_small_ehull_bilayers_mixed, fmt="%s")
    np.savetxt('monolayers_remaining_and_mixed_bilayers.npy', monolayer_remaining_bilayers_mixed, fmt="%s")

    
