"""
Unit and regression test for the WatCon package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

from WatCon.generate_static_networks import WaterNetwork
import MDAnalysis as mda
from MDAnalysis.topology import guessers
import numpy as np

def create_water_universe(n_residues: float) -> mda.Universe:
    """Create dummy water molecules to create WaterNetwork objects"""
    n_atoms = n_residues * 3
    resindices = np.repeat(range(n_residues), 3)
    segindices = [0] * n_residues

    # create empty universe
    sol = mda.Universe.empty(n_atoms, 
                             n_residues=n_residues,
                             atom_resindex=resindices, 
                             residue_segindex=segindices,
                             trajectory=True)

    atomindices = np.arange(1,n_atoms+1)
    sol.add_TopologyAttr('id', atomindices)
    sol.add_TopologyAttr('name', ['O', 'H1', 'H2']*n_residues)
    sol.add_TopologyAttr('type', ['O', 'H', 'H']*n_residues)
    sol.add_TopologyAttr('resname', ['SOL']*n_residues)
    sol.add_TopologyAttr('resid', list(range(1, n_residues+1)))
    sol.add_TopologyAttr('segid', ['SOL'])

    h2o = np.array([[ 0,        0,       0      ],  # oxygen
                     [ 0.95908, -0.02691, 0.03231],  # hydrogen
                     [-0.28004, -0.58767, 0.70556]]) # hydrogen

    grid_size = 5
    spacing = 3
    
    coordinates = []
    
    # translating h2o coordinates around a grid
    mol_ids = []
    for i in range(n_residues):
        x = spacing * (i % grid_size)
        y = spacing * ((i // grid_size) % grid_size)
        z = spacing * (i // (grid_size * grid_size))
    
        xyz = np.array([x, y, z])
    
        coordinates.extend(h2o + xyz.T)
        mol_ids.extend([i+1]*3)

    coords = np.array(coordinates)
    sol.atoms.positions = coords

    # add bonds
    bonds = []
    for o in range(0, n_atoms, 3):
        bonds.extend([(o, o+1), (o, o+2)])

    sol.add_TopologyAttr('bonds', bonds)
    sol.add_TopologyAttr('resids', mol_ids)
    return sol


class TestWaterNetwork:
    """Test WaterNetwork objects"""
    
    def test_add_water(self):
        u = create_water_universe(n_residues=4)
        wat_resid = u.select_atoms('resid 1').residues[0]
        net = WaterNetwork()
        oxygen = [atom for atom in wat_resid.atoms if 'O' in atom.name]
        hydrogens = [atom for atom in wat_resid.atoms if 'O' not in atom.name]

        net.add_water(wat_resid.resid, oxygen[0], wat_resid.resid, *hydrogens)
        assert len(oxygen)+len(hydrogens) == 3
        assert len(net.water_molecules) == 1

    def test_find_connections(self):
        u = create_water_universe(n_residues=4)
        # select only oxygens
        ag = u.select_atoms("all")
        net = WaterNetwork()
        #select only oxygen atoms
        for mol in ag.molecules:
            ats = [atom for atom in mol.atoms if 'O' in atom.name]
            net.add_water(mol.resid, *ats, mol.resid)

        connections = net.find_connections(dist_cutoff=3)
        assert len(connections) == 4
        connections = net.find_connectons(dist_cutoff=1)
        assert len(connections) == 0

        


        
#class TestWaterBox:
#    """Test simple pdb with eight waters"""
#    def __init__:
#        self.waterbox = "inputs/WATERS.pdb"
#        #Put other args in here
#
#    def get_connections(self.waterbox):
#        metrics, networks, _, _ = 
#
#    def check_connections(self.waterbox)
#
#    
#def test_water_connections(water_input):

