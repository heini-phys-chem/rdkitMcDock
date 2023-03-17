#!/usr/bin/env python3

import sys, os
import random
import warnings

import numpy as np
import openbabel
from openbabel import openbabel as ob

from rdkit import Chem

from utils_rdkit import *

# Suppress RuntimeWarning from first Metropolis-Hastings MC criterian, accept
warnings.filterwarnings("ignore")
os.system("rm -f out.xyz min.xyz directory/*.xyz")

def main():
    print("\nUsing OpenBabel version: {}\n------------------------------\n".format(openbabel.__version__))
    # Read cmd line arguments
    f_target, f_ligand, force_field, trajectories, steps, temperature = get_options(sys.argv[1:])

    # Read in molecules
    # Host = target
    target = readfile(f_target)
    # Guest = ligand
    ligand = readfile(f_ligand)
    print(" [+] Read molecules")

    # center ligand molecule
    com = get_com(ligand)
    print(" [+] get center of mass (ligand)")

    ligand = move_molecule(ligand, -com)
    print(" [+] move molecule to center (ligand)")

    # center target molecule
    com = get_com(target)
    print(" [+] get center of mass (target)")

    target = move_molecule(target, -com)
    print(" [+] move molecule to center (target) \n")

    numTot = target.GetNumAtoms() + ligand.GetNumAtoms()

    eb = calc_energy_xTB(ligand)
    ea = calc_energy_xTB(target)

    print()
    print("Energy Host before opt: {}".format(ea))
    print("Energy Guest before opt: {}".format(eb))
    print()

    ea, target = minimize_molecule_xTB(target)
    eb, ligand = minimize_molecule_xTB(ligand)

    print()
    print("Energy Host after opt: {}".format(ea))
    print("Energy Guest after opt: {}".format(eb))
    print()

    #writefile_xyz(ligand, "ligand_xtbopt_from_mcDock.xyz")
    #writefile(target, "target_xtbopt_from_mcDock.sdf")

    numTot = target.GetNumAtoms() + ligand.GetNumAtoms()
    e_low  = np.Inf
    eb_min = np.Inf


    c = 1
    numConfs = 1
    # Perfomr conformer search
#    print(" Perfomring rotor search for ligand molecule      file: {}".format(f_ligand))
#    ligand = set_conformations(ligand, force_field)
#
#    numConfs = ligand.NumConformers()
#
#    for c in range(1, numConfs+1):
#        ligand.SetConformer(c)
#
#        eb, ligand = minimize_molecule(ligand, force_field)
#
#        if eb < eb_min:
#            eb_min = eb
#
#        print(" Rotamere {}     E = {:.4f} kcal/mol".format(str(c), eb))
#
#    print(" Lowers energy conformation E = {:.4f} kcal/mol".format(eb_min))
#    print(" Running {} trajectories for {} steps".format(trajectories, steps))
#    print(" MC temperature (tau) = {}".format(temperature))


    # Start the MC simulation
    print("\n Conformation:         Trajectory:         Acceptance rate:               Final Ebind:")
    print(" -------------------------------------------------------------------------------------")

    # Loop over conformers
#    for c in range(1, numConfs+1):
#    ligand.SetConformer(c)

#    eb, ligand = minimize_molecule_xTB(ligand)

    # For every conformer loop over number of trajectories
    for n in range(int(trajectories)):
        # translate molecule
        direction = random_vector()
        com = get_com(ligand)
        temp = direction * 4.0 - com
        move_molecule(ligand, temp)

        # rotate molecule (Guest)
        rot = random_vector()
        theta = random_angle()# * np.pi
        rotate_molecule(ligand, rot, theta)

        # roll back to rejected MC moves
        # add both molecule to the same "system"
        mol_ligand = plus_equal_mols(target, ligand, numTot)
        mol_old    = mol_ligand
        #mol_old.SetCoordinates(mol_ligand.GetCoordinates())

#        ff = AllChem.MMFFGetMoleculeForceField(mol_ligand)
#        e = AllChem.UFFGetMoleculeForceField(mol_ligand).CalcEnergy()
        e = calc_energy_xTB(mol_ligand)

        energy_old = e
        delta_e = 0.0
        accept = 0

        startid = mol_ligand.GetNumAtoms() - ligand.GetNumAtoms()
        endid   = mol_ligand.GetNumAtoms()

        # MC simulation for every trajectory
        for step in range(int(steps)):
            move = random_vector()
            move *= random_length()
            move_molecule(mol_ligand, move, startid=startid, endid=endid)

            rot   = random_vector()
            theta = random_angle()
            rotate_molecule(mol_ligand, rot, theta, startid=startid, endid=endid)


            #e = AllChem.UFFGetMoleculeForceField(mol_ligand).CalcEnergy()

            e = calc_energy_xTB(mol_ligand)

            delta_e = e - energy_old

            # Metropolis-Hastings MC criterian, accept ...
            threash_hold = uniform1()
            if np.exp(-delta_e / float(temperature)) >= threash_hold:# and range_test:
                mol_old = mol_ligand
                energy_old = e
                accept += 1
            else:
                mol_ligand = mol_old
                e = energy_old


        ec, mol_ligand = minimize_molecule_xTB(mol_ligand)

        e_bind = ec - (ea + eb)

        acceptance_ratio = accept * 100.0 / (float(steps) + 1)
        acceptance_ratio = round(acceptance_ratio, 2)
        print(" {} / {}\t\t\t{} / {}\t\t\t{}\t%\t\t{:4.2f}\tkcal/mol".format(c, numConfs, str(n+1).zfill(2), trajectories, str(acceptance_ratio).zfill(2), e_bind), end = '')


        writefile_xyz(mol_ligand, "directory/out_{}_{}_{}.xyz".format(c, n, 1))


        if (e_bind < e_low):
            print("     <----- New Lowest")
            e_low = e_bind
            writefile_xyz(mol_ligand, "min.xyz")
        else:
            print()

    os.system("cat directory/out_*.xyz >> out.xyz")
    #os.system("rm -rf out_*.xyz")

    print("\n\n [+] Optimization done")

if __name__ == '__main__':
    main()
