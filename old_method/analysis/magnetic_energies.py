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
        energydiff = data["eDIFF"] 
        if 'energydiff' in locals():
            energydiff = energydiff/n

        name_energy.append(("/".join((str(folder).split("/")[-2:])), abs(energydiff) ))
        energies.append(energydiff)

    return energies, name_energy

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
        if os.path.isfile(f'{folder}/structure.json'):
        #if os.path.isfile(f'{folder}/results-asr.bilayersummary.json'):
            Ediff, name_energy = get_energy_difference(folder)
            magmoms, magmom = magnetic_type(folder)
            gap, name_material = get_gap_type(folder)

            #if not len(Ediff) == 0:
            #    if 0.05 > abs(Ediff[0]):
            #        if gap == 0.0:

            #            print(folder)
            #            print(Ediff[0])

            if not len(magmoms) == 0:
                exfoliation_energy, exfoliation_energy_material = get_exfoliation_energy(folder)
                ex_energy.extend(exfoliation_energy)
                
                print(folder)
                if gap == 0.0:
                    metals_exchange.append(Ediff)
                    metals_exfoliation.append(exfoliation_energy)
                else:
                    insulators_exchange.append(Ediff)
                    insulators_exfoliation.append(exfoliation_energy)

            Ediffs.extend(Ediff)
            name_energies.extend(name_energy)
  
    os.chdir('/home/niflheim2/cmr/WIP/stacking/tree-mads/')

    #np.savetxt('CrI3_AB_prime_insulators_exfoliation_energy_old.npy', insulators_exfoliation)   
    #np.savetxt('CrI3_AB_prime_insulators_exchange_energy_old.npy', insulators_exchange)
    #np.savetxt('CrI3_AB_prime_metals_exfoliation_energy_old.npy', metals_exfoliation)   
    #np.savetxt('CrI3_AB_prime_metals_exchange_energy_old.npy', metals_exchange)
 

    print('Number of insulators:', len(insulators_exchange))
    print('Number of metals:', len(metals_exchange))

    #plt.figure()
    #plt.scatter(metals_exfoliation, metals_exchange)
    #plt.scatter(insulators_exfoliation, insulators_exchange)
    #plt.xlabel('Exfoliation energy [eV]')
    #plt.ylabel('Energy difference pr. magnetic atom [eV/atom]')
    #plt.title('Energy_difference vs. Exfoliation energy')
    #plt.tight_layout()

    #plt.legend(['metals' ,'insulators'])
    #plt.savefig('Energy_diff_exfoliation_energy_scatterplot_colorcoded.png')

    #plt.figure()
    #plt.scatter(metals_exfoliation, metals_exchange)
    #plt.scatter(insulators_exfoliation, insulators_exchange)
    #plt.ylim([-0.0025, 0.0025])
    #plt.xlim([0, 0.04])
    #plt.xlabel('Exfoliation energy [eV]')
    #plt.ylabel('Energy difference pr. magnetic atom [eV/atom]')
    #plt.title('Energy_difference vs. Exfoliation energy')
    #plt.tight_layout()

    #plt.legend(['metals' ,'insulators'])
    #plt.savefig('CrI3_Energy_diff_exfoliation_energy_scatterplot_colorcoded_closeup.png')
        

    #plt.figure()
    #plt.scatter(metals_exfoliation, metals_exchange)
    #plt.scatter(insulators_exfoliation, insulators_exchange)
    #plt.ylim([-0.0001, 0.0001])
    #plt.xlim([0, 0.04])
    #plt.xlabel('Exfoliation energy [eV]')
    #plt.ylabel('Energy difference pr. magnetic atom [eV/atom]')
    #plt.title('Energy_difference vs. Exfoliation energy')
    #plt.tight_layout()

    #plt.legend(['metals' ,'insulators'])
    #plt.savefig('Energy_diff_exfoliation_energy_scatterplot_colorcoded_closer_closeup.png')

    #plt.figure()
    #plt.scatter(metals_exfoliation, metals_exchange, c='blue')
    #plt.ylim([-0.0001, 0.0001])
    #plt.xlim([0, 0.04])
    #plt.xlabel('Exfoliation energy [eV]')
    #plt.ylabel('Energy difference pr. magnetic atom [eV/atom]')
    #plt.title('Energy_difference vs. Exfoliation energy')
    #plt.tight_layout()

    #plt.legend(['metals'])
    #plt.savefig('CrI3_Energy_diff_exfoliation_energy_scatterplot_metals.png')

    #plt.figure()
    #plt.scatter(insulators_exfoliation, insulators_exchange, c='darkorange')
    #plt.ylim([-0.0001, 0.0001])
    #plt.xlim([0, 0.04])
    #plt.xlabel('Exfoliation energy [eV]')
    #plt.ylabel('Energy difference pr. magnetic atom [eV/atom]')
    #plt.title('Energy_difference vs. Exfoliation energy')
    #plt.tight_layout()

    #plt.legend(['insulators'])
    #plt.savefig('CrI3_Energy_diff_exfoliation_energy_scatterplot_insulators.png')

    #plt.figure()
    #plt.scatter(ex_energy,Ediffs)
    #plt.xlabel('Exfoliation energy')
    #plt.ylabel('Energy difference pr. magnetic atom [meV/atom]')
    #plt.title('Energy_difference vs. Exfoliation energy')
    #plt.tight_layout()

    #plt.savefig('Energy_diff_exfoliation_energy_scatterplot.png')


    #plt.figure()
    #plt.scatter(ex_energy,Ediffs)
    #plt.ylim([-0.2, 0.2])
    #plt.xlim([0, 0.04])
    #plt.xlabel('Exfoliation energy')
    #plt.ylabel('Energy difference pr. magnetic atom [meV/atom]')
    #plt.title('Energy_difference vs. Exfoliation energy')
    #plt.tight_layout()

    #plt.savefig('Energy_diff_exfoliation_energy_scatterplot_closeup.png')

