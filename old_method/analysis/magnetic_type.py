from pathlib import Path
from utils import listdirs
import os, re
from asr.core import read_json, ASRResult
import matplotlib.pyplot as plt
import numpy as np
from ase.io import read
from gpaw import GPAW


""" Run this scripts as:
python3 bilayerworkflow/magnetic_type.py tree/*/*/*/*/ """

def magnetic_type(folder):
    if os.path.exists(f"{folder}/results-asr.magstate.json"):
        data = read_json(f'{folder}/results-asr.magstate.json')
        magstate = data["magstate"]
        return magstate
    else:
        return "Missing"

def get_monolayer_data(monolayer_folder: str, fname: str = "results-asr.gs.json") -> ASRResult:
    path = "/".join([x for x in monolayer_folder.split("/") if x != ""][-3:])
    c2db_path = "/home/niflheim2/cmr/C2DB-ASR/tree/" + path + "/" + fname
    #assert Path(c2db_path).is_file(), c2db_path
    if  Path(c2db_path).is_file():
        data = read_json(c2db_path)
        magstate = str(data["magstate"])
        return magstate #read_json(c2db_path)
    else:
        return 'Missing'

def get_monolayer_id(monolayer_folder):
    monolayer_id = read_json(f"{monolayer_folder}/structure.json")
    monolayer_id_data = monolayer_id[1]
    data_id = monolayer_id_data["unique_id"]
    return data_id

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


    unique_monolayers = []

    total_number = []
    no_change = []
    NM_to_FM = []
    NM_to_AFM = []
    FM_to_AFM = []
    FM_to_NM = []
    AFM_to_FM = []
    AFM_to_NM = []

    for folder in folders:
        if os.path.isfile(f'{folder}/fm_U3.0.txt') and os.path.isfile(f'{folder}/afm_U3.0.txt') and os.path.isfile(f'{folder}/gs_U3.0.gpw'):
            
            deviation_matrix_fm, deviation_matrix_afm = check_magmoms(folder)

            monolayer_folder = re.split('[/]',f'{folder}')
            monolayer_folder.remove(monolayer_folder[-1])
            monolayer_folder = "/".join(monolayer_folder)

            magstate_bilayer = magnetic_type(folder)
            magstate_monolayer = get_monolayer_data(f"{monolayer_folder}", "results-asr.magstate.json")
        
            monolayer = get_monolayer_id(monolayer_folder)


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
            #    if AFM_norm == True:
                #exfoliation_energy, exfoliation_energy_material = get_exfoliation_energy(folder)
                #ex_energy.extend(exfoliation_energy)

                #if gap == 0.0:
                #    metals_exchange.append(Ediff)
                #    metals_exfoliation.append(exfoliation_energy)
                #else:
                #    insulators_exchange.append(Ediff)
                #    insulators_exfoliation.append(exfoliation_energy)

                if os.path.exists(f'{folder}/results-asr.bilayer_magnetism.json'):
                    total_number.append(folder)
        

                if f'{magstate_monolayer}' == 'NM' and f'{magstate_bilayer}' == 'FM':
                    NM_to_FM.append(folder)
                    print('NM to FM:', folder)
                if f'{magstate_monolayer}' == 'NM' and f'{magstate_bilayer}' == 'AFM':
                    NM_to_AFM.append(folder)
                    print('NM to AFM:', folder)
                if f'{magstate_monolayer}' == 'FM' and f'{magstate_bilayer}' == 'AFM':
                    FM_to_AFM.append(folder)
                if f'{magstate_monolayer}' == 'FM' and f'{magstate_bilayer}' == 'NM':
                    FM_to_NM.append(folder)
                if f'{magstate_monolayer}' == 'AFM' and f'{magstate_bilayer}' == 'FM':
                    AFM_to_FM.append(folder)
                if f'{magstate_monolayer}' == 'AFM' and f'{magstate_bilayer}' == 'NM':
                    AFM_to_NM.append(folder)
   
                if os.path.exists(f'{folder}/results-asr.bilayer_magnetism.json'):
                    if monolayer not in unique_monolayers:
                        unique_monolayers.append(monolayer)


    print('Total number of monolayers', len(unique_monolayers))
    print('Total number of bilayers:', len(total_number))
    print('No change:', len(no_change))
    print('NM to FM:', len(NM_to_FM))
    print('NM to AFM:', len(NM_to_AFM))
    print('FM to AFM:', len(FM_to_AFM))
    print('FM to NM:', len(FM_to_NM))
    print('AFM to FM:', len(AFM_to_FM))
    print('AFM to NM:', len(AFM_to_NM))
