from ase.io import read, write
from asr.core import read_json, command, option, DictStr, ASRResult
from asr.utils.bilayerutils import translation
import numpy as np
from ase.calculators.dftd3 import DFTD3
import os
from ase.visualize import view

def convert_mixer(mixer):
    from gpaw import MixerSum, MixerDif

    if 'type' not in mixer:
        raise ValueError('Key "type" is required')

    if mixer['type'] == 'MixerSum':
        beta = mixer['beta']
        nmaxold = mixer['nmaxold']
        weight = mixer['weight']
        return MixerSum(beta=beta, nmaxold=nmaxold, weight=weight)
    elif mixer['type'] == 'MixerDif':
        beta = mixer['beta']
        nmaxold = mixer['nmaxold']
        weight = mixer['weight']
        return MixerDif(beta=beta, nmaxold=nmaxold, weight=weight)
    else:
        raise ValueError(f'Unrecognized mixer type: {mixer["type"]}')


def check_magmoms(atoms, magmoms_FM, magmoms_AFM, magmom_AFM):
    
    name_material = []
    atoms = atoms
    
    magnetic_atoms = []
    magnetic_atom_types = []
    #TM3d_atoms = ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn']
    #TM3d_atoms = ['V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu']
    mag_elements = ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
                    'Y', 'Zr', 'Nb', 'Mo', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In',
                    'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl']

    for atom in atoms:
        if atom.symbol in mag_elements:
            magnetic_atoms.append(1)
            magnetic_atom_types.append(atom.symbol)
        else:
            magnetic_atoms.append(0)
            magnetic_atom_types.append("")
                
    #mag_atoms = [x for x, z in enumerate(magnetic_atoms) if z == 1]
    mag_atoms = []
    mag_atoms_type = []
    for iatom, z in enumerate(magnetic_atoms):
        if z==1: 
           mag_atoms.append(iatom)
           mag_atoms_type.append(magnetic_atom_types[iatom])

 
    magmoms_FM =  magmoms_FM
    magmoms_AFM = magmoms_AFM
    magmom_AFM = magmom_AFM

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
    
    ###############################
    check_values_fm = []
    
    #for x in deviation_matrix_fm:
    #    for y in x:
    #        if abs(y) > 0.05:
    #            check_values_fm.append(y)

    for m, x in zip(deviation_matrix_fm, mag_atoms_type):
        for n, y in zip(m, mag_atoms_type):
            if abs(n) > 0.1 and x == y:
                check_values_fm.append(n)

    if not len(check_values_fm) == 0:
        FM_state = 'ildefined'
    else:
        FM_state = 'FM'

    ##############################
    check_values_afm = []

    #if not deviation_matrix_afm == 'ildefined':
    #    for x in deviation_matrix_afm:
    #        for y in x:
    #            if abs(y) > 0.05:
    #                check_values_afm.append(y)

    if not deviation_matrix_afm == 'ildefined':
       for m, x in zip(deviation_matrix_afm, mag_atoms_type):
           for n, y in zip(m, mag_atoms_type):
               if abs(n) > 0.1 and x == y:
                   check_values_afm.append(n)


    if np.allclose((magmom_AFM/len(mag_atoms)),0, atol=0.01):
        deviation_matrix_afm == 'ildefined'

    if deviation_matrix_afm == 'ildefined':
        check_values_afm.append('ildefined')

    if not len(check_values_afm) == 0:
        AFM_state = 'ildefined'
    else:
        AFM_state = 'AFM'

    return FM_state, AFM_state



def get_bilayer(atoms, top_layer, magmoms, d=None, config=None):
    assert config is not None, 'Must provide config'

    # Extract data necessary to construct the bilayer from
    # the separate layers. This allows us to easily set
    # the magnetic moments in the separate layers.
    translation_data = read_json('translation.json')
    if 'Saddle_from' in translation_data:
       z1 = top_layer.positions[0,2]
       saddle_atoms = read('structure.json')
       z2 = saddle_atoms.positions[len(top_layer),2]
       h = abs(z2-z1)
       print('optimal height for saddle',h, z1, z2)
    elif d is not None and d>=0:
       h = d
    else:
       if os.path.isfile("results-asr.zscan.json") or os.path.isfile('zscan_correction.txt'):
          relax_data = read_json("results-asr.zscan.json")
          h = relax_data["optimal_height"]
       elif os.path.isfile("results-asr.relax_bilayer.json"):
          relax_data = read_json("results-asr.relax_bilayer.json")
          h = relax_data["optimal_height"]

    t_c = np.array(read_json('translation.json')[
                   'translation_vector']).astype(float)

    # Make local copies so we dont modify input
    atoms = atoms.copy()
    top = top_layer.copy()

    # Set the magnetic moments
    if config.lower() == 'fm':
        atoms.set_initial_magnetic_moments(magmoms)
        top.set_initial_magnetic_moments(magmoms)
    elif config.lower() == 'afm':
        atoms.set_initial_magnetic_moments(magmoms)
        top.set_initial_magnetic_moments(-magmoms)
    else:
        raise ValueError(f'Configuration not recognized: {config}')

    # Combine the two layers
    bilayer = translation(t_c[0], t_c[1], h, top, atoms)

    return bilayer

#@command(module='asr.interlayer_magnetic_exchange',
#         requires=['results-asr.magstate.json'],
#         resources='24:10h',
#         dependencies=['asr.magstate',
#                       'asr.relax_bilayer'])

#@command(module='asr.interlayer_magnetic_exchange',
#         resources='24:10h',
#         dependencies=['asr.relax_bilayer'])

@command(module='asr.interlayer_magnetic_exchange',
         resources='24:10h')
@option('-c', '--calculator', help='Calculator params.', type=DictStr())
@option('-u', '--hubbardu', type=float, default=None, help="Hubbard U correction")  # U = 3
@option('-d', '--interlayer', type=float, default=None, help="Interlayer distance")
@option('-ml', '--mlfolder', type=str, help="monolayer folder")
#@option('-m', '--mixer', help='help', type=DictStr())
def main(calculator: dict = {
        'name': 'gpaw',
        'mode': {'name': 'pw', 'ecut': 800},
        'xc': 'PBE',
        'basis': 'dzp',
        'kpts': {'density': 12.0, 'gamma': True},
        'occupations': {'name': 'fermi-dirac',
                        'width': 0.05},
        'maxiter': 5000,
        'mixer': {"method": "sum",
        #'mixer': {'type': 'mixerdif',  # used for second try
                  "beta": 0.02,
                  "history": 5,
                  "weight": 50},
        'poissonsolver' : {'dipolelayer': 'xy'},
        'convergence': {'bands': 'CBM+3.0', "energy": 1e-6, "density": 1e-6},
        'nbands': '200%'},
        #mixer: dict = None,
        interlayer: float = -100.0, #if negative read from zscan
        mlfolder: str = '..',
        hubbardu: float = 0) -> ASRResult:

    """Calculate the energy difference between FM and AFM configurations.

    Returns the energy difference between the FM and
    AFM configurturations for a bilayer, where FM is
    defined as all spins pointing in the same directions,
    and AFM is defined as spins in the two layers
    pointing in opposite directions.

    If U is non-zero, Hubbard-U correction is used for 3d TM
    atoms, if there are any.

    If mixer is not None, a custom mixer is used for the GPAW calculation.
    """


    #if mixer is not None:
    #    raise NotImplementedError  # "params" is not defined
    #    params['mixer'] = convert_mixer(mixer)

    u = hubbardu
    atoms = read(f'{mlfolder}/structure.json')
    top_layer = read('toplayer.json')

    # Fix the cell of the bottom and top layers
    # to be the bilayer cell
    cell = read(f"{mlfolder}/structure.json").cell
   
    maxz = np.max(atoms.positions[:, 2])
    minz = np.min(atoms.positions[:, 2])
    w = maxz - minz
    vacuum = 6
    cell[2, 2] += vacuum + w

    atoms.cell = cell
    top_layer.cell = cell

    # Read the magnetic moments from the monolayers.
    # This is used as the starting point for the bilayer
    # magnetic moments.
    magmoms = read_json(f"{mlfolder}/structure.json")[1]["magmoms"]

    if hubbardu is None:
        # 3d TM atoms which need a Hubbard U correction
        TM3d_atoms = {'V':3.1, 'Cr':3.5, 'Mn':3.8, 'Fe':4.0, 'Co':3.3, 'Ni':6.4, 'Cu':4.0}
        atom_ucorr = set([atom.symbol for atom in atoms if atom.symbol in TM3d_atoms])
        U_corrections_dct = {symbol: f':d, {TM3d_atoms[symbol]}' for symbol in atom_ucorr}
    if hubbardu is not None:
        u = hubbardu
        TM3d_atoms = ['V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu']
        atom_ucorr = set([atom.symbol for atom in atoms
                      if atom.symbol in TM3d_atoms])
        U_corrections_dct = {symbol: f':d, {u}' for symbol in atom_ucorr}


    mag_elements = ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
                    'Y', 'Zr', 'Nb', 'Mo', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In',
                    'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl']

    magnetic_atoms = []
    for atom in atoms:
        if atom.symbol in mag_elements:
            magnetic_atoms.append(1)
        else:
            magnetic_atoms.append(0)
                
    mag_atoms = [x for x, z in enumerate(magnetic_atoms) if z == 1]
    for x in mag_atoms:
        magmoms[x] = abs(magmoms[x])

    calculator.update(setups=U_corrections_dct)

    from ase.calculators.calculator import get_calculator_class
    name = calculator.pop('name')
    calc_fm = get_calculator_class(name)(**calculator, txt=f"fm_U.txt") #had {u} before
    calc_afm = get_calculator_class(name)(**calculator, txt=f"afm_U.txt") #had {u} before

    # We use cutoff=60 for the vdW correction to be consistent with
    calc_fm_D3 = DFTD3(dft=calc_fm, cutoff=60)

    # FM Calculation
    bilayer_FM = get_bilayer(atoms, top_layer, magmoms, d=interlayer, config="FM")
    view(bilayer_FM)
    initial_magmoms_FM = bilayer_FM.get_initial_magnetic_moments()

    bilayer_FM.set_calculator(calc_fm)
    final_magmoms_FM = bilayer_FM.get_magnetic_moments()

    bilayer_FM.set_calculator(calc_fm_D3)
    eFM = bilayer_FM.get_potential_energy()

    FM_syms = bilayer_FM.get_chemical_symbols()
    for x in [i for i, e in enumerate(FM_syms) if e in mag_elements]:
        #if abs(final_magmoms_FM[x])>0.1 or abs(initial_magmoms_FM[x])>0.1:
        #   assert np.sign(final_magmoms_FM[x]) == np.sign(initial_magmoms_FM[x])
 
        mag_max = max(abs(np.array(initial_magmoms_FM))) 
        if abs(initial_magmoms_FM[x])>(0.1*mag_max):
           assert np.sign(final_magmoms_FM[x]) == np.sign(initial_magmoms_FM[x])

    # AFM Calculation
    calc_afm_D3 = DFTD3(dft=calc_afm, cutoff=60)

    bilayer_AFM = get_bilayer(atoms, top_layer, magmoms, d=interlayer, config="AFM")
    bilayer_AFM.set_calculator(calc_afm_D3)

    eAFM = bilayer_AFM.get_potential_energy()

    eDIFF = eFM - eAFM

    magmoms_fm = calc_fm.get_magnetic_moments()
    M_FM = calc_fm.get_magnetic_moment()

    magmoms_afm = calc_afm.get_magnetic_moments()
    M_AFM = calc_afm.get_magnetic_moment()

    monolayer_dipole = read_json(f'{mlfolder}/results-asr.gs.json')['dipz']
    if abs(monolayer_dipole)<1e-3:
       print('Magmoms checked not to change much')
       FM_state, AFM_state = check_magmoms(atoms, magmoms_fm, magmoms_afm, M_AFM)
    else: 
       print('Magmoms not compared due to out-of-plane dipole')
       FM_state = 'FM'
       AFM_state = 'AFM'

    ### Keep the structure before reinitialising the magnetic moments
    #atoms0 = read(f'structure.json')
    #atoms0.write(f'structure_beforemag.json')

    if FM_state == 'FM' and AFM_state == 'AFM':
        if eDIFF > 0:
            calc_afm.write(f'gs_U.gpw') #had {u} before
            bilayer_AFM.write(f'structure.json')
        else:
            calc_fm.write(f'gs_U.gpw') #had {u} before
            bilayer_FM.write(f'structure.json')

    if FM_state == 'FM' and AFM_state == 'ildefined':
        calc_fm.write(f'gs_U.gpw')  #had {u} before
        bilayer_FM.write(f'structure.json')
    
    if FM_state == 'ildefined' and AFM_state == 'AFM':
        calc_afm.write(f'gs_U.gpw')  #had {u} before
        bilayer_AFM.write(f'structure.json')


    return dict(eFM=eFM, eAFM=eAFM, eDIFF=eDIFF,
                M_AFM=M_AFM, M_FM=M_FM,
                magmoms_fm=magmoms_fm, magmoms_afm=magmoms_afm,
                FM_state=FM_state, AFM_state=AFM_state)

if __name__ == '__main__':
    main.cli()
