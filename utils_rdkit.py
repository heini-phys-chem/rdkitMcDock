import sys, os, getopt

import numpy as np
from openbabel import openbabel as ob

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D

from ase.io import read, write
from ase.optimize import BFGS
from xtb.ase.calculator import XTB

import xyz2mol as x2m


def get_options(argv):
   '''
   Read in cmd line arguments
   ligand: ligand geometry (sdf)
   target: target geometry (sdf)
   forceField: force field to be used
   trajectories: how many trajectorues (MC runs) per conformer
   steps: number of MC steps
   temperature: temperature for acceptance test
   '''
   #opts, args = getopt.getopt(argv, "hl:t:f:j:s:u:m:", ["help", "ligand=", "target=", "forceField=", "trajectories=", "steps=", "temperature=", "mutations="])
   opts, args = getopt.getopt(argv, "hl:t:f:j:s:u:", ["help", "ligand=", "target=", "forceField=", "trajectories=", "steps=", "temperature="])
   target, ligand, ff, trajectories, steps, temperature = "PA_hexamer.sdf", "Ru_guest_for_PA.sdf", "UFF", 5, 200, 0.3
   for opt, arg in opts:
      if opt == '--help':
         print("How to run:")
         print('test.py --ligand <ligand> --target <target> --forceField <force field> --trajectories # --steps # --temperature\n')
         print(" -> ligand and target molecules as .sdf files")
         print(" -> force fileds as implemented in openbabel (UFF, MMFF94, ...)")
         print(" -> trajectories:\t# of trajectories")
         print(" -> steps:\tnumber of MC steps")
         print(" -> temperature for acceptance step\n")
         print(" -> mutations: # of mutations\n")
         sys.exit()
      elif opt in ("-l", "--ligand"):
         ligand = arg
      elif opt in ("-t", "--target"):
         target = arg
      elif opt in ("-f", "--forceField"):
         ff = arg
      elif opt in ("-j", "--trajectories"):
         trajectories = arg
      elif opt in ("-s", "--steps"):
         steps = arg
      elif opt in ("-u", "--temperature"):
         temperature = arg
#      elif opt in ("-m", "--mutations"):
#          mutations = arg

   #return target, ligand, ff, trajectories, steps, temperature, mutations
   return target, ligand, ff, trajectories, steps, temperature


def readfile(f):
    '''
    input:   filename (sdf)
    returns: rdkit mol object
    '''
    mol = Chem.MolFromMolFile(f, removeHs=False, sanitize=True)
    mol = Chem.RWMol(mol)
    return mol


def writefile(mol, f):
    '''
    input:   mol object, filename (sdf)
    returns: -
    '''
    print(Chem.MolToMolBlock(mol), file=open(f, 'w+'))

def writefile_xyz(mol, f):
    numAtoms = mol.GetNumAtoms()
    labels = []
    coords = []

    NAME = {1:"H", 6:"C", 7:"N", 8:"O", 44:"Ru"}
    for i, atom in enumerate(mol.GetAtoms()):
        label = atom.GetAtomicNum()
        labels.append(label)
        positions = mol.GetConformer().GetAtomPosition(i)
        coords.append([positions.x, positions.y, positions.z])

    labels = [NAME[label] for label in labels]

    fout = open(f, 'w')
    fout.write("{}\n\n".format(numAtoms))

    for i, label in enumerate(labels):
        fout.write("{} {} {} {}\n".format(label, coords[i][0], coords[i][1], coords[i][2]))


    fout.close()

def get_com(mol):
    '''
    get center of mass
    '''
    com = np.array([0.0, 0.0, 0.0])

    for i in range(mol.GetNumAtoms()):
        positions = mol.GetConformer().GetAtomPosition(i)
        tmp = np.array([ positions.x, positions.y, positions.z ])
        com += tmp

    com /= mol.GetNumAtoms()

    return com

def move_molecule(mol, move, startid=0, endid=-1):
    '''
    move molecule using move vector
    '''

    if endid == -1:
        endid = mol.GetNumAtoms()

    for i in range(startid, endid):
        positions = mol.GetConformer().GetAtomPosition(i)
        temp = np.array([positions.x, positions.y, positions.z])
        temp += move
        mol.GetConformer().SetAtomPosition(i, Point3D(temp[0], temp[1], temp[2]))

    return mol

def minimize_molecule(mol, force_field):
    '''
    optimizes molecule using a force field
    returns optimized geometry and energy
    '''
    res = AllChem.UFFOptimizeMoleculeConfs(mol)
    e = res[0][1]

    return e, mol

def minimize_molecule_xTB(mol):
    cmd = "ulimit -s unlimited"
    os.system(cmd)
    writefile_xyz(mol, "run_opt.xyz")
    mol_ase = read('run_opt.xyz')

    calc = XTB(method="GFN2-xTB")
    mol_ase.calc = calc

    dyn = BFGS(mol_ase)
    dyn.run(fmax=0.05)

    e = mol_ase.get_potential_energy()

    write("run_min.xyz", mol_ase)

    cmd = "python ~/workcopies/xyz2mol/xyz2mol.py run_min.xyz -o sdf > run_min.sdf"
    os.system(cmd)

    mol_min = readfile("run_min.sdf")

    cmd = "rm -f run_opt.xyz run_min.xyz run_min.sdf"
    os.system(cmd)

    return e, mol_min

def calc_energy_xTB(mol):
    cmd = "ulimit -s unlimited"
    os.system(cmd)
    writefile_xyz(mol, "run_sp.xyz")
    mol_ase = read('run_sp.xyz')

    calc = XTB(method="GFN2-xTB")
    mol_ase.calc = calc

    e = mol_ase.get_potential_energy()

    cmd = "rm -f xtbrestart xtbtopo.mol run_sp.xyz wbo charges"
    os.system(cmd)

    return e



def rotate(V, J, T):
    x = V[0]
    y = V[1]
    z = V[2]

    u = J[0]
    v = J[1]
    w = J[2]

    norm = np.sqrt(u*u + v*v + w*w)
    inv_norm_sqrt = 1.0 / (norm*norm)

    sint = np.sin(T)
    cost = np.cos(T)

    a = (u * (u*x + v*y + w*z) + (x * (v*v + w*w) - u * (v*y + w*z)) * cost + norm * (-w*y + v*z) * sint) * inv_norm_sqrt
    b = (v * (u*x + v*y + w*z) + (y * (u*u + w*w) - v * (u*x + w*z)) * cost + norm * ( w*x - u*z) * sint) * inv_norm_sqrt
    c = (w * (u*x + v*y + w*z) + (z * (u*u + v*v) - w * (u*x + v*y)) * cost + norm * (-v*x + u*y) * sint) * inv_norm_sqrt

    rotated = np.array([a, b, c])

    return rotated

def rotate_molecule(mol, direction, theta, startid=0, endid=-1):
    '''
    rotate molecule using rot matrix, theta
    '''
    if endid == -1:
        endid = mol.GetNumAtoms()

    com = np.array([0.0, 0.0, 0.0])

    for i in range(startid, endid):
        positions = mol.GetConformer().GetAtomPosition(i)
        com += np.array([positions.x, positions.y, positions.z])

    com /= (endid - startid)

    for i in range(startid, endid):
        positions = mol.GetConformer().GetAtomPosition(i)
        temp = np.array([positions.x, positions.y, positions.z])
        temp -= com
        mol.GetConformer().SetAtomPosition(i, Point3D(temp[0], temp[1], temp[2]))
        temp = rotate(temp, direction, theta)
        temp += com
        mol.GetConformer().SetAtomPosition(i, Point3D(temp[0], temp[1], temp[2]))

def random_vector():
    '''
    generates random 3D vector [-1, 1]
    '''
    rng = np.random.default_rng()
    vec = rng.random((3,))
    x = np.random.choice([-1., 1.], size=1)
    y = np.random.choice([-1., 1.], size=1)
    z = np.random.choice([-1., 1.], size=1)

    vec_tmp = np.array([ x*vec[0], y*vec[1], z*vec[2] ])
    vec = vec_tmp.reshape(3,)

    return vec

def uniform1():
    '''
    returns random number [0, 1] (scalar)
    '''
    rng = np.random.default_rng()
    return rng.random()

def random_angle():
    '''
    returns random angle for molecule rotation (scalar)
    '''
    rng = np.random.default_rng()
    return 90./180.*np.pi * rng.random()

def random_length():
    '''
    returns random length [0, 0.3] (scalar)
    '''
    rng = np.random.default_rng()
    #return 2 * rng.random()
    ranInt = rng.integers(low=0, high=30, size=1)


    if ranInt == 0:
        return 0
    else:
        return ranInt / 100.

def set_conformations(mol, force_field):
    '''
    Conformer search
    '''
    ff = ob.OBForceField.FindForceField(force_field)
    ff.Setup(mol)

    rmsd_cutoff   = 0.5
    energy_cutoff = 50.0
    conf_cutoff   = 10000000
    verbose       = False

    ff.DiverseConfGen(rmsd_cutoff, conf_cutoff, energy_cutoff, verbose)
    ff.GetConformers(mol)


    return mol

def plus_equal_mols(target, ligand, numTot):
    '''
    mol += mol2 does not exist in the python wrapper from openbabel.
    Therefore, this function adds the ligand to the target molecule (same mol object).
    In the first step (# of atoms of target < # num atoms target+ligand) it is simply added
    by creating the ligand atoms and bonds.
    For every subsequent step, the ligand atom positions are simply updated
    '''
    #ligand = ob.OBMol()
    #ligand = readfile_xyz("Ru_guest_for_PA.xyz")

    if target.GetNumAtoms() < numTot:
        mol_new = target
        numAtoms = target.GetNumAtoms()

        atom0      = []
        atom1      = []
        bond_order = []

        for bond in ligand.GetBonds():
            atom0.append(bond.GetBeginAtomIdx()+numAtoms)
            atom1.append(bond.GetEndAtomIdx()+numAtoms)
            bond_order.append(bond.GetBondType())

        for i, atom in enumerate(ligand.GetAtoms()):
            atomic_num = atom.GetAtomicNum()
            positions = ligand.GetConformer().GetAtomPosition(i)
            temp = np.array([positions.x, positions.y, positions.z])

            mol_new.AddAtom(Chem.Atom(atomic_num))
            mol_new.GetConformer().SetAtomPosition(i+numAtoms, Point3D(temp[0], temp[1], temp[2]))

        for i in range(len(atom0)):
            mol_new.AddBond(atom0[i], atom1[i], bond_order[i])

    else:
        mol_new = target
        numAtoms = target.GetNumAtoms() - ligand.GetNumAtoms()

        for i in range(ligand.GetNumAtoms()):
            positions = ligand.GetConformer().GetAtomPosition(i)
            X = positions.x
            Y = positions.y
            Z = positions.z

            mol_new.GetConformer().SetAtomPosition(numAtoms+i, Point3D(X, Y, Z))

    Chem.SanitizeMol(mol_new)

    return mol_new

