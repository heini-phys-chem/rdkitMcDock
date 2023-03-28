# rdkitMcDock

This code takes a `Host` molecule and a `Guest` molecule, reads them in, optimizes the geometries, merge them into the same `system`, moves both to the origin, and then performs
an MC simulation by rotating the `Guest` molecule and optimizing (L-BFGS) both geometries using ASE and the xTB calculator.

## Installation:
Conda needs to be installed on your system. Then:

- create a conda environment (here it's called McDock)
- activate the environment
- install requirements (openbabel, RdKit, ASE, xTB (with python wrappers), numpy, xyz2mol)
- xyz2mol also needs to be present locally (cloned from github) to fix a problem with reading xyz files

copy paste the commands after the $ sign to your terminal

```
(base)   $ conda create -name McDock
(base)   $ conda activate McDock
(McDock) $ conda install -c anaconda numpy
(McDock) $ conda install -c conda-forge openbabel
(McDock) $ conda install -c conda-forge xtb xtb-python
(McDock) $ conda install -c conda-forge ase
(McDock) $ conda install -c conda-forge rdkit
(McDock) $ conda install -c conda-forge xyz2mol
```

To fix and use xyz2mol clone the git repository:
```
git clone https://github.com/jensengroup/xyz2mol
```
cd into the directory
```
cd xyz2mol
```
Now, change line 550 of xyz2mol.py from
```
550                 atomic_symbol, x, y, z = line.split()
```
to
```
550                 atomic_symbol, x, y, z = line.split()[:4]
```
This is to fix a problem where the energies are included in the xyz file, causing the split to fail. We can safely ignore the energies since we recompute them anyway.

Finally, type `make`
```
make
```

The xyz2mol path is hardcoded wherever you cloned it. To run McDock you need to change line 158 in utils_rdkit.py to the path where you cloned the xyz2mol repository.
```
158     cmd = "python ~/workcopies/xyz2mol/xyz2mol.py run_min.xyz -o sdf > run_min.sdf"
```
change `~/workcopies/` to the path where you cloned xyz2mol!

## Run mcDock
Before running you need to change the memory default options of xTB using following command in your cmd line and also make a directory called directory:
```
ulimit -s unlimited
mkdir directory
```
Then you can run mcDock using following command:
```
./mcDock_rdkit.py --ligand small_conformer.sdf --target PA_hexamer.sdf  --forceField UFF --trajectories 3 --steps 15 --temperature 0.3
```
forcefield is obsolete but at the moment necessary...

## mcDock Workflow
![alt text](figs/mcDock-workflow.png)
