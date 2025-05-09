"""
Generate water networks based on dynamical data
"""

import os, sys
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import distances
from joblib import Parallel, delayed  # For parallel processing
import networkx as nx
import matplotlib.pyplot as plt

from scipy.spatial import cKDTree
from scipy.spatial.distance import pdist, squareform
import numpy as np

from WatCon.sequence_processing import *
import WatCon.sequence_processing as sequence_processing
from WatCon.visualize_structures import project_clusters
import WatCon.residue_analysis as residue_analysis

class WaterAtom:
    """
    Object for storing information of atoms in water molecules

    Attributes
    ----------
    index : int
        PDB atom index
    coordinates : tuple
        Coordinates of atom
    resname : str
        Name of residue. Default is WAT.
    name: str
        Name of atom
    resid: int
        Residue number of atom
    """
    def __init__(self, index, atom_name, residue_number, x, y, z):
        """
        Initialize the class

        Parameters
        ----------
        index : int
            PDB atom index
        atom_name: str
            Name of atom
        residue_number: int
            Residue number of atom
        x: float
            x-coordinate of position
        y: float
            y-coordinate of position
        z: float
            z-coordinate of position
        """
        self.index = index
        self.coordinates = (x, y, z)
        self.resname = 'WAT'
        self.name = atom_name
        self.resid = residue_number

class WaterMolecule:
    """
    Object for storing information of water molecules

    Attributes
    ----------
    index : int
        PDB atom index
    H1 : WaterAtom
        WaterAtom object for first hydrogen
    H2 : WaterAtom
        WaterAtom object for second hydrogen
    O : WaterAtom
        WaterAtom object for oxygen
    resname : str
        Name of residue. Default is WAT.
    resid: int
        Residue number of water molecule
    """
    def __init__(self, index, O: WaterAtom, H1: WaterAtom, H2: WaterAtom, residue_number):
        """
        Initialize the class

        Parameters
        ----------
        index : int
            PDB atom index
        O: WaterAtom
            WaterAtom object for oxygen atom
        H1: WaterAtom 
            WaterAtom object for first hydrogen
        H2: WaterAtom
            WaterAtom obejct for second hydrogen
        residue_number: int
            Residue number of atom
        """
        self.index = index
        self.H1 = H1
        self.H2 = H2
        self.O = O
        self.resname = 'WAT'
        self.resid = residue_number

class OtherAtom:
    """
    Object for storing information of other atoms (not waters)

    Attributes
    ----------
    index : int
        PDB atom index
    coordinates : tuple
        Coordinates of atom
    resname : str
        Name of residue.
    msa_resid: int
        Common residue index from MSA    
    resid: int
        Residue number of atom    
    name: str
        Name of atom
    """
    def __init__(self, index, atom_name, residue_name, x, y, z, residue_number, msa_residue_number, hbonding):
        """
        Initialize the class

        Parameters
        ----------
        index : int
            PDB atom index
        atom_name: str
            Name of atom
        residue_name: str
            Name of residue
        x: float
            x-coordinate of position
        y: float
            y-coordinate of position
        z: float
            z-coordinate of position
        residue_number: int
            Residue number of atom
        """
        self.index = index
        self.coordinates = (x, y, z)
        self.resname = residue_name
        self.msa_resid = msa_residue_number
        self.resid = residue_number
        self.name = atom_name
        self.hbonding = hbonding
  
class WaterNetwork:  #For water-protein analysis -- extrapolate to other solvent maybe
    """
    Object for storing information regarding water-water and water-protein connections

    Attributes
    ----------
    water_molecules : list
        List of WaterMolecule objects
    protein_atoms : list
        List of OtherAtom objects
    connections : list
        List of connections among atoms
    active_region: list
        Collection of WaterMolecule and OtherAtom objects making up a user-defined active site   
    graph: NetworkX graph object
        NetworkX graph object containing nodes and edges as defined by self.connections 
    """
    def __init__(self):
        """
        Initialize the class
        """
        self.water_molecules = []
        self.protein_atoms = []
        #self.protein_subset = []
        #self.molecules = []
        self.connections = None
        self.active_region = None
        self.graph = None

    def add_atom(self, index, atom_name, residue_name, x, y, z, residue_number=None, msa_residue_number=None):
        """
        Add atom (OtherAtom object)

        Parameters
        ----------
        index : int
            PDB atom index
        atom_name: str
            Name of atom
        residue_name : str
            Name of residue.
        x : float
            x-coordinate of atom
        y : float
            y-coordinate of atom
        z : float
            z-coordinate of atom
        residue_number: int
            Residue number of atom  
        msa_residue_number: int
            Common residue index from MSA    

        Returns
        ----------
        None
        """
        if residue_name == 'WAT' or residue_name == 'HOH' or residue_name == 'SOL': #these names are hardcoded
            water  = WaterMolecule(index, atom_name, x, y, z) 
            self.water_molecules.append(water)
        else: #Currently only adds protein atoms
            if ('O' in atom_name) or ('N' in atom_name) or ('S' in atom_name) or ('P' in atom_name):
                hbonding = True
            else:
                hbonding = False

            mol = OtherAtom(index, atom_name, residue_name, x, y, z, residue_number, msa_residue_number, hbonding)
            self.protein_atoms.append(mol)

    def add_water(self, index, o, h1, h2, residue_number):
        """
        Add water molecule

        Parameters
        ----------
        index : int
            PDB atom index
        atom_name: str
            Name of atom
        o: MDAnalysis atom object
            MDAnalysis atom object corresponding to water oxygen
        residue_number: int
            Residue number of water  
        h1: MDAnalysis atom object
            MDAnalysis atom object corresponding to first water hydrogen
        h2: MDAnalysis atom object
            MDAnalysis atom object corresponding to second water hydrogen

        Returns
        ----------
        None
        """
        o = WaterAtom(o.index, 'O', residue_number, *o.position)
        h1 = WaterAtom(h1.index, 'H1',residue_number, *h1.position)
        h2 = WaterAtom(h2.index, 'H2',residue_number, *h2.position)
        water = WaterMolecule(index, o, h1, h2, residue_number)
        self.water_molecules.append(water)

    def select_active_region(self, reference, box, active_region_radius=8.0, active_region_COM=False):
        """
        Select active site atoms based on distance to reference atoms.

        Parameters
        ----------
        reference : MDAnalysis AtomGroup object
            MDAnalysis Atom group (or list of atoms) defining the reference.
        box : array-like
            Simulation box used for periodic boundary conditions (PBC).
        dist_cutoff : float, optional
            Distance cutoff for selection. Default is 8.0 Å.
        active_region_COM : bool, optional  
            Whether to take center of mass of active site references, or combine for sphere selection
        Returns
        -------
        tuple
            - self.active_region : list  
              Selected active site atoms.  
            - active_region_protein : list  
              Active site protein atoms.  
            - active_region_water : list  
              Active site water atoms.  
        """

        #Create empty set for active site atoms
        active_region_atoms = []


        if active_region_COM is False:
            #Find coordinates for refrence point
            reference_positions = np.array([ref.position for ref in reference])  # Precompute reference positions

        else:
            reference_positions = reference.center_of_mass()


        #Find protein atoms in active site
        protein_active = []

        reference_resids = {r.resid for r in reference}  # Set of reference resids for fast lookup
        for atm in self.protein_atoms:
            #Immediately include atoms which are a part of the reference
            if atm.resid in reference_resids:
                active_region_atoms.append(atm)
                protein_active.append(atm)

            #Include atoms within a distance cutoff
            else:
                dist = np.min(distances.distance_array(np.array(atm.coordinates).reshape(1, -1), reference_positions, box=box))
                if dist <= active_region_radius:
                    active_region_atoms.append(atm)
                    protein_active.append(atm)

        #Find water molecules in active site
        water_active = []

        for mol in self.water_molecules:
            #Include atoms within a distance cutoff
            water_positions = np.array([mol.O.coordinates, mol.H1.coordinates, mol.H2.coordinates])
            
            dist = np.min(distances.distance_array(water_positions, reference_positions))
            if dist <= active_region_radius:            
                active_region_atoms.append(mol)
                water_active.append(mol)


        self.active_region = list(active_region_atoms)  # Convert set back to list if order matters
        return self.active_region, list(protein_active), list(water_active)

    def find_connections(self, dist_cutoff=3.3, water_active=None, protein_active=None, active_region_only=False, water_only=False):
        """
        Find the shortest connections using a K-D Tree.

        Parameters
        ----------
        dist_cutoff : float, optional
            Distance cutoff for connections. Default is 3.3 Å.
        water_active : list
            Selection of active site water molecules.
        protein_active : list
            Selection of active site protein molecules.
        active_region_only : bool, optional
            If True, only find connections among active site atoms. Default is False.
        water_only : bool, optional
            If True, only find connections among waters. Default is False.

        Returns
        -------
        list of tuples
            Each connection is represented as a tuple with the following elements:
            - connections[0] : int  
                Index of the first atom.
            - connections[1] : int  
                Index of the second atom.
            - connections[2] : str  
                Atom name of the first atom.
            - connections[3] : str  
                Type of interaction ('WAT-WAT' or 'WAT-PROT').
            - connections[4] : str  
                Whether the interaction is in the active site ('active_region') or not ('not_active_region').
            - connections[5] : str
                Protein atom type ('backbone' or 'sidechain')
        """

        connections = []

        # Select active site atoms if specified
        if active_region_only:
            waters = water_active
            protein = protein_active
        else:
            waters = self.water_molecules
            protein = self.protein_atoms

        # Extract water oxygen coordinates and indices
        water_coords = np.array([mol.O.coordinates for mol in waters])
        water_indices = np.array([mol.O.index for mol in waters])
        water_names = np.array(['O' for _ in waters])
        

        # Water-Water connections
        tree = cKDTree(water_coords)
        dist, indices = tree.query(water_coords, k=10, distance_upper_bound=dist_cutoff)

        for i, neighbors in enumerate(indices):
            for j, neighbor in enumerate(neighbors):

                if neighbor != i and dist[i, j] <= dist_cutoff:
                    site_status = (
                        'None' if self.active_region is None else
                        'active_region' if any(waters[i].resid == f.resid for f in self.active_region) else
                        'not_active_region'
                    )
                    if water_indices[i] < water_indices[neighbor]:
                        connections.append((water_indices[i], water_indices[neighbor], water_names[i], 'WAT-WAT', site_status))

        # Water-Protein connections
        if not water_only:
            protein_coords = np.array([atm.coordinates for atm in protein])
            protein_indices = np.array([atm.index for atm in protein])
            protein_names = np.array([atm.name for atm in protein])

            tree = cKDTree(protein_coords)
            dist, indices = tree.query(water_coords, k=10, distance_upper_bound=dist_cutoff)

            for i, neighbors in enumerate(indices):
                for j, neighbor in enumerate(neighbors):
                    if dist[i, j] <= dist_cutoff:
                        site_status = (
                            'None' if self.active_region is None else
                            'active_region' if any(waters[i].resid == f.resid or protein[neighbor].resid == f.resid for f in self.active_region) else
                            'not_active_region'
                        )
                        if protein_names[neighbor] == 'O' or protein_names[neighbor] == 'N':
                            classification = 'backbone'
                        else:
                            classification = 'side-chain'
                        connections.append((protein_indices[neighbor], water_indices[i], protein_names[neighbor], 'WAT-PROT', site_status, classification))

        return connections


    def find_directed_connections(self, dist_cutoff=2.5, water_active=None, protein_active=None, active_region_only=False, 
                                    water_only=False, angle_criteria=None):
        """
        Find the shortest directed connections using a K-D Tree.

        Parameters
        ----------
        dist_cutoff : float, optional
            Distance cutoff for connections. Default is 2.0 Å.
        water_active : list
            Selection of active site water molecules.
        protein_active : list
            Selection of active site protein molecules.
        active_region_only : bool, optional
            If True, only find connections among active site atoms. Default is False.
        water_only : bool, optional
            If True, only find connections among waters. Default is False.
        angle_criteria: float, optional
            Additional angle criteria for defining hydrogen bonds. Default is None.

        Returns
        -------
        list of tuples
            Each connection is represented as a tuple with the following elements:
            - connections[0] : int  
              Index of the first atom.
            - connections[1] : int  
              Index of the second atom.
            - connections[2] : str  
              Atom name of the first atom.
            - connections[3] : str  
              Type of interaction ('WAT-WAT' or 'WAT-PROT').
            - connections[4] : str  
              Whether the interaction is in the active site ('active_region') or not ('not_active_region').
        """

        # Initialize empty list for connections
        connections = []

        # Select active site atoms if specified
        if active_region_only:
            waters = water_active
            protein = protein_active
        else:
            waters = self.water_molecules
            protein = self.protein_atoms

        # Gather water indices and coordinates
        water_H_indices = []
        water_H_coords = []
        water_H_names = []

        water_O_coords = []
        water_O_indices = []
        water_O_names = []

        for mol in waters:

            #Select status
            if self.active_region is None:
                site_status = 'None'
            elif mol in self.active_region:
                site_status = 'active_region'
            else:
                site_status = 'not_active_region'

            # Add H1 atom
            water_H_indices.append(mol.O.index)  #Use only O index
            water_H_coords.append(mol.H1.coordinates)
            water_H_names.append('H1')

            # Add H2 atom
            water_H_indices.append(mol.O.index) #Use only O index
            water_H_coords.append(mol.H2.coordinates)
            water_H_names.append('H2')

            # Add O
            water_O_indices.append(mol.O.index)
            water_O_coords.append(mol.O.coordinates)
            water_O_names.append('O')

        # Convert water_coords to a numpy array
        water_H_coords = np.array(water_H_coords).reshape(-1,3)
        water_O_coords = np.array(water_O_coords).reshape(-1,3)

        #Find protein-water connections
        if water_only == False:

            #Gather indices and coords
            protO_coords = []
            protO_indices = []        
            protO_names = []

            protH_coords = []
            protH_indices = []
            protH_names = []

            for atm in protein:
                if 'P' in atm.name or 'O' in atm.name or 'N' in atm.name or 'S' in atm.name:
                    protO_coords.append(atm.coordinates)
                    protO_indices.append(atm.index)
                    protO_names.append(atm.name)

                elif 'H' in atm.name:
                    protH_coords.append(atm.coordinates)
                    protH_indices.append(atm.index)
                    protH_names.append(atm.name)

            protH_coords = np.array(protH_coords).reshape(-1,3)
            protO_coords = np.array(protO_coords).reshape(-1,3)

            #Find distances between protein H and water O
            #Create KDTree using protein-H coordinates
            tree = cKDTree(protH_coords)

            #Query for distances with water O coordinates
            dist, indices = tree.query(water_O_coords, k=10, distance_upper_bound=dist_cutoff)

            for index_near, index_ref in enumerate(indices):
                for i, distance in enumerate(dist[index_near]):
                    #Check for cutoff, scipy will output infs if distance is too high
                    if distance <= dist_cutoff:
                        if not active_region_only:
                            # Determine active site status
                            if self.active_region is None:
                                site_status = 'None'
                            elif any(idx in [water_H_indices[index_near], protH_indices[index_ref[i]]] for idx in [f.index for f in self.active_region]):
                                site_status = 'active_region'
                            else:
                                site_status = 'not_active_region'
                        else: 
                            site_status = 'active_region'

                        if protH_names[index_ref[i]] == 'H' or protH_names[index_ref[i]] == 'HA': #MAKE THIS BETTER
                            classification = 'backbone'
                        else:
                            classification = 'side-chain'
                        #Append connections
                        if angle_criteria is None:
                            connections.append([protH_indices[index_ref[i]], water_O_indices[index_near], protH_names[index_ref[i]] , 'WAT-PROT', site_status, classification])
                        else:
                            protein_hydrogen_coords = protH_coords[index_ref[i]]
                            water_oxygen_coords = water_O_coords[index_near]

                            protein_resid = [atm.resid for atm in self.protein_atoms if atm.index == protH_indices[index_ref[i]]][0]

                            protein_O_coordinates = [atm.coordinates for atm in self.protein_atoms if (atm.resid == protein_resid and 'H' not in atm.name)]
                            distances = [np.linalg.norm(protein_hydrogen_coords-f) for f in protein_O_coordinates]
                            arg = np.argmin(distances)

                            prot_heavy_coordinates = protein_O_coordinates[arg]

                            prot_water = water_oxygen_coords - protein_hydrogen_coords
                            prot_prot = prot_heavy_coordinates - protein_hydrogen_coords

                            prot_water = prot_water.flatten()
                            prot_prot = prot_prot.flatten()

                            cosine_angle = np.dot(prot_water, prot_prot) / (np.linalg.norm(prot_water) * np.linalg.norm(prot_prot))
                            angle1 = np.degrees(np.arccos(cosine_angle))

                            if angle1 >= angle_criteria:
                                connections.append([protH_indices[index_ref[i]], water_O_indices[index_near], protH_names[index_ref[i]] , 'WAT-PROT', site_status, classification])

            #Find distances between protein O,S,P,N and water H

            #Create KDTree using protein OSPN coordinates
            tree = cKDTree(protO_coords)

            #Query for distances with water-H coords
            dist, indices = tree.query(water_H_coords, k=10, distance_upper_bound=dist_cutoff)
            for index_near, index_ref in enumerate(indices):
                for i, distance in enumerate(dist[index_near]):
                    if distance <= dist_cutoff:        
                        if not active_region_only:
                            # Determine active site status
                            if self.active_region is None:
                                site_status = 'None'
                            elif any(idx in [water_H_indices[index_near], protO_indices[index_ref[i]]] for idx in [f.index for f in self.active_region]):
                                site_status = 'active_region'
                            else:
                                site_status = 'not_active_region'
                        else: 
                            site_status = 'active_region'

                        if protO_names[index_ref[i]] == 'O' or protO_names[index_ref[i]] == 'N':
                            classification = 'backbone'
                        else:
                            classification = 'side-chain'
                        #Append connections
                        if angle_criteria is None:
                            connections.append([water_H_indices[index_near], protO_indices[index_ref[i]], water_H_names[index_near], 'WAT-PROT', site_status,classification])
                        else:
                            protein_heavy_coords = protO_coords[index_ref[i]]
                            water_hydrogen_coords = water_H_coords[index_near]
                            water_o_coords = [water.O.coordinates for water in self.water_molecules if (water.O.index == water_H_indices[index_near])]

                            prot_water = protein_heavy_coords - water_o_coords
                            water1 = water_o_coords - water_hydrogen_coords

                            prot_water = prot_water.flatten()
                            water1 = water1.flatten()

                            cosine_angle = np.dot(prot_water, water1) / (np.linalg.norm(prot_water) * np.linalg.norm(water1))
                            angle1 = np.degrees(np.arccos(cosine_angle))

                            if angle1 >= angle_criteria:
                                connections.append([water_H_indices[index_near], protO_indices[index_ref[i]], water_H_names[index_near], 'WAT-PROT', site_status, classification])


            
        #Find distances between water O and water H

        #Create KDTree for water-O coords
        tree = cKDTree(water_O_coords)

        #Query for distances with water-H coords
        dist, indices = tree.query(water_H_coords, k=10, distance_upper_bound=dist_cutoff)

        
        for index_near, index_ref in enumerate(indices):
            for i, distance in enumerate(dist[index_near]):

                if distance <= dist_cutoff:
                    if not active_region_only:
                        # Determine active site status
                        if self.active_region is None:
                            site_status = 'None'
                        elif any(idx in [water_H_indices[index_near], water_O_indices[index_ref[i]]] for idx in [f.index for f in self.active_region]):
                            site_status = 'active_region'
                        else:
                            site_status = 'not_active_region'
                    else: 
                        site_status = 'active_region'


                    if water_H_indices[index_near] != water_O_indices[index_ref[i]]: #Check to make sure connection is not within the same water
                        
                        #Append connections
                        if angle_criteria is None:

                            if water_H_indices[index_near]<water_O_indices[index_ref[i]]:
                                #print('IMPORTANT: CHECKING DUPLICATE CONNECTIONS')
                                connections.append([water_H_indices[index_near],water_O_indices[index_ref[i]], water_H_names[index_near], 'WAT-WAT', site_status])

                        else:
                            water_hydrogen_coords = water_H_coords[index_near]
                            water_o1_coords = water_O_coords[index_ref[i]]
                            
                            water_o2_coords = [water.O.coordinates for water in self.water_molecules if (water.O.index == water_H_indices[index_near])][0]
                            water1 = water_hydrogen_coords - water_o1_coords
                            water2 = water_hydrogen_coords - water_o2_coords

                            water1 = water1.flatten()
                            water2 = water2.flatten()

                            cosine_angle = np.dot(water1, water2) / (np.linalg.norm(water2) * np.linalg.norm(water1))
                            angle1 = np.degrees(np.arccos(cosine_angle))

                            if angle1 >= angle_criteria and water_H_indices[index_near]<water_O_indices[index_ref[i]]:
                                connections.append([water_H_indices[index_near],water_O_indices[index_ref[i]], water_H_names[index_near], 'WAT-WAT', site_status])


        #Find distances between water H and water O

        #Create KDTree for water-O coords
        tree = cKDTree(water_H_coords)

        #Query for distances with water-H coords
        dist, indices = tree.query(water_O_coords, k=10, distance_upper_bound=dist_cutoff)

        
        for index_near, index_ref in enumerate(indices):
            for i, distance in enumerate(dist[index_near]):

                if distance <= dist_cutoff:
                    if not active_region_only:
                        # Determine active site status
                        if self.active_region is None:
                            site_status = 'None'
                        elif any(idx in [water_O_indices[index_near], water_H_indices[index_ref[i]]] for idx in [f.index for f in self.active_region]):
                            site_status = 'active_region'
                        else:
                            site_status = 'not_active_region'
                    else: 
                        site_status = 'active_region'

                    if water_O_indices[index_near] != water_H_indices[index_ref[i]]: #Check to make sure connection is not within the same water
                        
                        #Append connections
                        if angle_criteria is None:

                            if water_O_indices[index_near]<water_H_indices[index_ref[i]]:
                                #print('IMPORTANT: CHECKING DUPLICATE CONNECTIONS')
                                if [water_H_indices[index_ref[i]], water_O_indices[index_near], water_H_names[index_ref[i]], 'WAT-WAT', site_status] not in connections:
                                    connections.append([water_H_indices[index_ref[i]], water_O_indices[index_near], water_H_names[index_ref[i]], 'WAT-WAT', site_status])

                        else:
                            water_hydrogen_coords = water_H_coords[index_ref[i]]
                            water_o1_coords = water_O_coords[index_near]
                            
                            water_o2_coords = [water.O.coordinates for water in self.water_molecules if (water.O.index == water_H_indices[index_ref[i]])][0]
                            water1 = water_hydrogen_coords - water_o1_coords
                            water2 = water_hydrogen_coords - water_o2_coords

                            water1 = water1.flatten()
                            water2 = water2.flatten()

                            cosine_angle = np.dot(water1, water2) / (np.linalg.norm(water2) * np.linalg.norm(water1))
                            angle1 = np.degrees(np.arccos(cosine_angle))

                            #if angle1 >= angle_criteria and water_H_indices[index_ref[i]]<water_O_indices[index_near]:
                            if angle1 >= angle_criteria and [water_H_indices[index_ref[i]], water_O_indices[index_near], water_H_names[index_ref[i]], 'WAT-WAT', site_status] not in connections and [water_O_indices[index_near], water_H_names[index_ref[i]], 'WAT-WAT', site_status] not in connections:
                                connections.append([water_H_indices[index_ref[i]], water_O_indices[index_near], water_H_names[index_ref[i]], 'WAT-WAT', site_status])

                

        return connections


    def generate_oxygen_network(self, box, msa_indexing=None, active_region_reference=None, active_region_COM=False, active_region_radius=8.0, 
                                active_region_only=False, water_only=False, max_connection_distance=3.0):
        """
        Generate network based only on oxygens -- direct comparability to static structure networks

        Parameters
        ----------
        msa_indexing : list, optional
            list of MSA residue indexes, default is None
        active_region_reference: MDAnalysis AtomGroup object, optional
            Reference to define active site. Default is None.
        active_region_COM : bool, optional  
            Whether to take center of mass of active site references, or combine for sphere selection
        active_region_radius: float, optional
            Radius to define active site from active_region_reference. Default is 8.0 Å
        active_region_only : bool, optional
            If True, only find connections among active site atoms. Default is False.
        water_only : bool, optional
            If True, only find connections among waters. Default is False.
        max_connection_distance: float, optional
            Distance cutoff for connections. Default is 2.0 Å.


        Returns
        -------
        self.graph
            NetworkX graph object describing the network.
        """

        #Initialize graph
        G = nx.Graph()

        #Select active site
        if active_region_reference is not None:
            self.active_region, protein_active, water_active = self.select_active_region(active_region_reference, box=box, active_region_radius=active_region_radius, active_region_COM=active_region_COM)


        #Use MSA indexing
        if msa_indexing is not None:
            MSA_indices = msa_indexing
        else:
            MSA_indices = ['X'] * 10000 #Dummy list


        #Only include atoms in active site -- greatly increases performance
        if active_region_only:
            for molecule in water_active:
                G.add_node(molecule.O.index, pos=molecule.O.coordinates, atom_category='WAT', MSA=None) #have nodes on all oxygens

            if water_only == False:
                for molecule in protein_active:   
                    MSA_index = MSA_indices[molecule.resid-1]           
                    G.add_node(molecule.index, pos=molecule.coordinates, atom_category='PROTEIN', MSA=MSA_index)

            self.connections = self.find_connections(dist_cutoff=max_connection_distance, water_active=water_active, protein_active=protein_active, active_region_only=active_region_only, water_only=water_only)
            for connection in [f for f in self.connections if f[4]=='active_region']:
                G.add_edge(connection[0], connection[1], connection_type=connection[3], active_region=connection[4])

        #Include all atoms
        else:
            for molecule in self.water_molecules:
                G.add_node(molecule.O.index, pos=molecule.O.coordinates, atom_category='WAT', MSA=None) #have nodes on all oxygens

            if water_only == False:
                #for molecule in self.protein_subset:
                for molecule in self.protein_atoms:
                    MSA_index = MSA_indices[molecule.resid-1]  
                    G.add_node(molecule.index, pos=molecule.coordinates, atom_category='PROTEIN', MSA=None)
            
            self.connections = self.find_connections(dist_cutoff=max_connection_distance, water_active=None, protein_active=None, active_region_only=False, water_only=water_only)

            for connection in self.connections:
                G.add_edge(connection[0], connection[1], connection_type=connection[3], active_region=connection[4])

        #Save as self.graph
        self.graph = G

        return self.graph

    def generate_directed_network(self, box, msa_indexing=None, active_region_reference=None, active_region_COM=False, active_region_radius=8.0, 
                                  active_region_only=False, water_only=False, angle_criteria=None, max_connection_distance=2.0):
        """
        Generate directed graph using H -> O directionality

        Parameters
        ----------
        msa_indexing : list, optional
            list of MSA residue indexes, default is None
        active_region_reference: MDAnalysis AtomGroup object, optional
            Reference to define active site. Default is None.
        active_region_COM : bool, optional  
            Whether to take center of mass of active site references, or combine for sphere selection
        active_region_radius: float, optional
            Radius to define active site from active_region_reference. Default is 8.0 Å
        active_region_only : bool, optional
            If True, only find connections among active site atoms. Default is False.
        water_only : bool, optional
            If True, only find connections among waters. Default is False.
        angle_criteria: float, optional
            Additional angle criteria to define hydrogen bonds. Defualt is None
        max_connection_distance: float, optional
            Distance cutoff for connections. Default is 2.0 Å.


        Returns
        -------
        self.graph
            NetworkX graph object describing the network.
        """

        G = nx.DiGraph() 
 
        #Select active site if a reference is given
        if active_region_reference is not None:
            self.active_region, protein_active, water_active = self.select_active_region(active_region_reference, box=box, active_region_radius=active_region_radius, active_region_COM=active_region_COM)

        #Use MSA indexing
        if msa_indexing is not None:
            MSA_indices = msa_indexing
        else:
            MSA_indices = ['X'] * 10000 #Dummy list

        #Only active site atoms in networks
        if active_region_only==True:
            #Add nodes
            for molecule in water_active:
                G.add_node(molecule.O.index, pos=molecule.O.coordinates, atom_category='WAT', MSA=None) #have nodes on all oxygens

            if water_only == False:
                for molecule in protein_active:              
                    MSA_index = MSA_indices[molecule.resid-1]
                    G.add_node(molecule.index, pos=molecule.coordinates, atom_category='PROTEIN', MSA=MSA_index)

            #Add edges
            self.connections = self.find_directed_connections(dist_cutoff=max_connection_distance, water_active=water_active, protein_active=protein_active, active_region_only=active_region_only, water_only=water_only, angle_criteria=angle_criteria)
            for connection in [f for f in self.connections if f[4]=='active_region']:
                G.add_edge(connection[0], connection[1], connection_type=connection[3], active_region=connection[4])

        #All atoms in network
        else:
            #Add nodes
            for molecule in self.water_molecules:
                G.add_node(molecule.O.index, pos=molecule.O.coordinates, atom_category='WAT', MSA=None) #have nodes on all oxygens

            if water_only == False:
                #for molecule in self.protein_subset:
                for molecule in self.protein_atoms:
                    G.add_node(molecule.index, pos=molecule.coordinates, atom_category='PROTEIN', MSA=MSA_index)
            
            #Add edges
            self.connections = self.find_directed_connections(dist_cutoff=max_connection_distance, water_active=None, protein_active=None, active_region_only=False, water_only=water_only)
            for connection in self.connections:
                G.add_edge(connection[0], connection[1], connection_type=connection[3], active_region=connection[4])

        self.graph = G
        return G

    def get_density(self, selection='all'):
        """
        Calculate the density of the graph.

        Requires `self.graph` to exist. The density is calculated as:

        .. math::
            \\text{density} = \\frac{N_{edges}}{(N_{nodes} * (N_{nodes} - 1) / 2)}

        This represents the ratio between actual edges and possible edges.

        Parameters
        ----------
        selection : {'all', 'active_region', 'not_active_region'}
            Specifies which subset of the graph to analyze.

        Returns
        -------
        float
            The density of the selected network.
        """   
        #Choose all subgraphs under particular criteria
        if selection=='all':
            S = self.graph
        else:
            S = self.graph.edge_subgraph([(edge1, edge2) for (edge1,edge2, data) in S.edges(data=True) if data['active_region']==selection])

        #Calculate density for subgraph
        nedges = S.number_of_edges()
        nnodes = S.number_of_nodes()
        density = nedges/((nnodes*(nnodes-1))/2) #density is ratio between edges and possible edges

        return density

    def get_connected_components(self, selection='all'):
        """
        Compute the connected components of the graph.

        Requires `self.graph` to exist. Uses `weakly_connected_components` if the graph is directed, 
        and `connected_components` if the graph is undirected.

        Parameters
        ----------
        selection : {'all', 'active_region', 'not_active_region'}
            Specifies which subset of the graph to analyze.

        Returns
        -------
        numpy.ndarray
            An array of connected components in the selected network.
        """
        #Initiate empty array for connected components
        components = []

        #Choose all subgraphs under particular criteria
        if selection=='all':
            S = self.graph
        else:
            S = self.graph.edge_subgraph([(edge1, edge2) for (edge1,edge2, data) in self.graph.edges(data=True) if data['active_region']==selection])

        #Use weakly_connected_components for directed graph
        if self.graph.is_directed():
            for val in [len(cc) for cc in nx.weakly_connected_components(S)]:
                components.append(val)
        #Use connected_components for undirected
        else:
            for val in [len(cc) for cc in nx.connected_components(S)]:
                components.append(val)

        #Reshape components to make plotting easier
        components = np.array(components).reshape(-1,1)
        return components
    
    def get_interactions(self, selection='all'):
        """
        Get a dictionary of water-water and protein-water interactions.

        Parameters
        ----------
        selection : {'all', 'active_region', 'not_active_region'}
            Specifies which subset of the graph to analyze.

        Returns
        -------
        dict
            A dictionary where keys are interaction types (e.g., 'WAT-WAT', 'WAT-PROT') 
            and values are the number of occurrences.
        """
        interaction_dict = residue_analysis.get_interaction_counts(self, selection)
        return(interaction_dict)
    
    def get_per_residue_interactions(self, selection='all'):
        """
        Get a dictionary of interactions per residue.

        Parameters
        ----------
        selection : {'all', 'active_region', 'not_active_region'}
            Specifies which subset of the graph to analyze.

        Returns
        -------
        dict
            A dictionary where keys are residue identifiers and values 
            are the number of interactions per residue.
        """
        residue_interaction_dict = residue_analysis.get_per_residue_interactions(self, selection)
        return(residue_interaction_dict)
    
    def get_CPL(self, selection='all', calculate_path='all', exclude_single_points=False):
        """
        Calculate the characteristic path length (CPL).

        The CPL is the average shortest path length in the network. 
        Components of length 1 (isolated nodes) are not included in the calculation.

        Parameters
        ----------
        selection : {'all', 'active_region', 'not_active_region'}
            Specifies which subset of the graph to analyze.
        exclude_single_points : bool, optional
            If True, excludes isolated nodes from the calculation. Default is False.

        Returns
        -------
        float
            The characteristic path length of the selected network.
        """
        #Choose all subgraphs under particular criteria
        if selection=='all':
            S = self.graph
        else:
            S = self.graph.edge_subgraph([(edge1, edge2) for (edge1,edge2, data) in self.graph.edges(data=True) if data['active_region']==selection])

        try:
            CPL = nx.average_shortest_path_length(S)
        except nx.NetworkXError: #average_shortest_path_length will fail if graph is not connected
            if S.is_directed():

                #Change to undirected graph for nx.average_shortest_path_length
                S = S.to_undirected() 

            if exclude_single_points:
                cc = [f for f in nx.connected_components(S) if len(f)>1]
            else:
                cc = [f for f in nx.connected_components(S)]

            CPLs = []
            if calculate_path == 'all':
                for C in (S.subgraph(c).copy() for c in cc):
                    CPLs.append(nx.average_shortest_path_length(C))

                #Average over all calculated CPLs
                CPL = np.array(CPLs).mean()
            
            else:
                largest_cc = max(nx.connected_components(S), key=len)
                CPL = nx.average_shortest_path_length(largest_cc)

        return CPL

    def get_shortest_path(self, selection='all', source=None, target=None):
        """
        Calculate the shortest path for a given network

        Parameters
        ----------
        selection : {'all', 'active_region', 'not_active_region'}
            Specifies which subset of the graph to analyze.
        source : int, optional
            Source node to initialize path
        target : int, optional
            Target node to terminate path

        Returns
        -------
        list
            Nodes corresponding to the shortest path.


        Note
        ----
        Node indexes are equivalent to MDAnalysis 0-based atom indexing (Add 1 to your pdb atom numbering!)
        """
        if selection=='all':
            S = self.graph
        else:
            S = self.graph.edge_subgraph([(edge1, edge2) for (edge1,edge2, data) in self.graph.edges(data=True) if data['active_region']==selection])

        shortest_path = nx.shortest_path(S, source, target)
        return shortest_path
    
    def get_clustering_coefficient(self, selection='all'):
        """
        Calculate the clustering coefficient for each node

        Parameters
        ----------
        selection : {'all', 'active_region', 'not_active_region'}
            Specifies which subset of the graph to analyze.

        Returns
        -------
        dict
            Clustering coefficient at each node position      
        """
        #Choose all subgraphs under particular criteria
        if selection=='all':
            S = self.graph
        else:
            S = self.graph.edge_subgraph([(edge1, edge2) for (edge1,edge2, data) in self.graph.edges(data=True) if data['active_region']==selection])

        CC_dict = nx.clustering(S)
        return CC_dict

    def get_entropy(self, selection='all'):
        """
        Calculate the graph entropy.

        The method is adapted from:
        https://stackoverflow.com/questions/70858169/networkx-entropy-of-subgraphs-generated-from-detected-communities

        Parameters
        ----------
        selection : {'all', 'active_region', 'not_active_region'}
            Specifies which subset of the graph to analyze.

        Returns
        -------
        float
            The computed graph entropy.
        """
        def degree_distribuiton(G):
            vk = dict(G.degree())
            vk = list(vk.values())

            maxk = np.max(vk)
            mink = np.min(vk)

            kvalues = np.arange(0, maxk+1)

            Pk = np.zeros(maxk+1)
            for k in vk:
                Pk[k] = Pk[k] + 1

            Pk = Pk/sum(Pk)
            return kvalues, Pk
        
        if selection=='all':
            S = self.graph
        else:
            S = self.graph.edge_subgraph([(edge1, edge2) for (edge1,edge2, data) in self.graph.edges(data=True) if data['active_region']==selection])

        k, Pk = degree_distribuiton(S)

        H = 0
        for p in Pk:
            if p > 0:
                H = H - p*np.log2(p)

        return H
    
    def get_all_coordinates(self, selection='all', water_only=True):
        """
        Retrieve all atomic coordinates, useful for clustering.

        Parameters
        ----------
        selection : {'all', 'active_region', 'not_active_region'}
            Specifies which subset of atoms to include.
        water_only : bool, optional
            If True, only includes water molecule coordinates. If False, includes both 
            water and interacting protein atoms. Default is False.

        Returns
        -------
        numpy.ndarray
            An array of coordinates for the selected atoms.
        """
        #Choose all subgraphs under particular criteria
        if selection=='all':
            #Find all coordinates -- only water oxygens
            coords = [np.array(f.O.coordinates) for f in self.water_molecules]

            if not water_only:
                print('Including OtherAtoms in clustering')
                #coords.extend([np.array(f.coordinates) for f in self.protein_subset])
                coords.extend([np.array(f.coordinates) for f in self.protein_atoms])

        else:
            #Find all coordinates -- only water oxygens
            coords = [np.array(f.O.coordinates) for f in self.active_region if type(f)==WaterMolecule]

            if not water_only:
                coords.extend([np.array(f.coordinates) for f in self.active_region if type(f)==OtherAtom])


        return coords


def get_clusters(coordinates, cluster, min_samples=10, eps=0.0, n_jobs=1, filename_base='DYNAMIC_CLUSTER'):
    """
    Get clusters over a series of WaterNetwork objects.

    Parameters
    ----------
    coordinates : array-like
        Array of coordinates combined for all networks
    cluster : {'hdbscan', 'dbscan', 'optics'}
        Clustering method
    min_samples : int, optional
        Minimum number of points per cluster. Default is 10
    coordinates : numpy.ndarray, optional
        Collection of coordinates
    eps: float, optional
        eps value for neighborhood size when clustering. Default is 0.0
    n_jobs: int, optional
        Available cores for parallelization. Default is 1 (no parallelization)
    filename_base: str, optional
        Filename base for PDB of cluster centers

    Returns
    ----------
    tuple
        - List of cluster labels
        - Dictionary describing cluster_centers
    """
    from WatCon.find_conserved_networks import cluster_coordinates_only
    from WatCon.visualize_structures import project_clusters
    cluster_labels, cluster_centers = cluster_coordinates_only(coordinates, cluster, min_samples, eps, n_jobs)
    project_clusters(cluster_centers, filename_base=filename_base)

    return cluster_labels, cluster_centers

def collect_densities(topology_file, trajectory_file, active_region_definition, active_region_COM, custom_selection, 
                      water_name, water_oxygen, output_name):
    """
    Calculate density of water positions, output density file, calculate hotspots from densities, and output PDB file corresponding to cluster hotspots.

    Parameters
    ----------
    pdb_file : str
        Full path to MDAnalysis-readable topology file
    trajectory_file : str
        Full path to MDAnalysis-readable trajectory file
    active_region_definition : str
        MDAnalysis selection language to define active site
    active_region_COM : bool
        Indicate whether to use center of mass for active site definition
    custom_selection : str
        MDAnalysis selection language to include custom residues in protein definition
    water_name : str
        Resname for water
    water_oxygen : str
        Atom name for water oxygens
    output_name : str
        Basename for output file (will automatically be given the .dx extension)
    
    Returns
    -------
    array-like
        Coordinates of density hotspots
    """
    from MDAnalysis.analysis import align, density
    from find_conserved_networks import find_clusters_from_densities

    #Perform density analysis
    ref = mda.Universe(topology_file)
    if not isinstance(trajectory_file, list):
        trajectory_file = [trajectory_file]
    u = mda.Universe(topology_file, *trajectory_file)
    align.AlignTraj(u, ref, select='name CA', filename='tmp.dcd', in_memory=False).run()

    u = mda.Universe(topology_file, 'tmp.dcd')
    ag = u.select_atoms(f'protein or {custom_selection}')
    ag.write(f'{output_name}.pdb')

    if active_region_COM:
        reference = u.select_atoms('resid 220 or resid 214').center_of_mass()
        x, y, z = reference  # Unpack coordinates
        print(x,y,z)
        # Properly format the selection string
        ow = u.select_atoms(f"(resname {water_name} and name {water_oxygen}) and point {x} {y} {z} 9", updating=True)
    
    else:
        ow = u.select_atoms(f"(resname {water_name} and name {water_oxygen}) and {active_region_definition}", updating=True)

    D = density.DensityAnalysis(ow, delta=1.0)
    D.run()
    D.results.density.convert_density('TIP3P')
    D.density.export(f"{output_name}.dx", type='double')

    hotspot_coords = find_clusters_from_densities(f"{output_name}.dx", output_name=f"{output_name}_densityclusters", threshold=1.5)
    return hotspot_coords

def extract_objects_per_frame(pdb_file, trajectory_file, frame_idx, network_type, custom_selection, 
                              active_region_reference, active_region_COM, active_region_radius, water_name, msa_indexing, 
                              active_region_only=False, directed=False, angle_criteria=None, max_connection_distance=3.0):
    """
    Extract and compute a water network for each frame.

    This function initializes a network based on the provided parameters and returns a 
    `WaterNetwork` object representing the computed network.

    Parameters
    ----------
    pdb_file : str
        Path to the topology file.
    trajectory_file: str
        Path to the trajectory file
    frame_idx: int
        Index of given frame
    network_type : {'water-water', 'water-protein'}
        Type of network to construct.
    custom_selection : str or None
        MDAnalysis selection string for custom residue selections.
    active_region_reference : str or None
        MDAnalysis selection string defining the reference for the active site.
    active_region_COM : bool, optional  
        Whether to take center of mass of active site references, or combine for sphere selection
    active_region_radius : float
        Radius (in Å) to define the active site region.
    water_name : str
        Name of water molecules in the system.
    msa_indexing : bool
        Whether to use MSA (multiple sequence alignment) indexing.
    active_region_only : bool, optional
        If True, only includes active site atoms in the network. Default is False.
    directed : bool, optional
        If True, constructs a directed network. Default is False.
    angle_criteria : float or None, optional
        Angle cutoff criteria for hydrogen-bonding structures. Default is None.
    max_connection_distance : float, optional
        Maximum distance (in Å) for defining connections in the network. Default is 3.0.

    Returns
    -------
    WaterNetwork
        A `WaterNetwork` object representing the computed network for the given PDB.
    """
 
    #Allow for custom residues in protein selection
    if custom_selection is None:
        custom_sel = ''
    else:
        custom_sel = f"or {custom_selection}" 

    #Allow for user-defined water name
    if water_name is None:
        water = 'resname HOH or resname WAT or resname SOL or resname H2O'
    else:
        water = f"resname {water_name}"

    #Create mda Universe
    u = mda.Universe(pdb_file, trajectory_file) 

    try:
        #Find maximum distance between edge of protein and middle of protein
        protein = u.select_atoms("protein")
        protein_com = protein.center_of_mass()
        distances_to_com = np.linalg.norm(protein.atoms.positions - protein_com, axis=1)
        max_distance = np.max(distances_to_com)

        #Separate key atom groups
        ag_wat = u.select_atoms(f'{water} and (sphzone {max_distance+0.5} protein)', updating=True)
        if not directed:
            ag_protein = u.select_atoms(f'(protein {custom_sel}) and (name N* or name O* or name P* or name S*)', updating=True)
        else:
            # Restrict hydrogens to those near polar atoms to minimize guess_bonds() overhead
            polar_heavy = u.select_atoms(f'(protein {custom_sel}) and (name N* or name O* or name P* or name S*)')
            hydrogens = u.select_atoms(f'protein {custom_sel} and name H*')

            # Find hydrogens near these heavy atoms (within 1.2 Å, a typical H-bond distance)
            from MDAnalysis.analysis.distances import distance_array

            dists = distance_array(hydrogens.positions, polar_heavy.positions)
            close_hydrogens = hydrogens[dists.min(axis=1) < 1.2]  # Select only close hydrogens

            # Step 4: Combine hydrogens and polar atoms into one AtomGroup
            relevant_atoms = close_hydrogens + polar_heavy
            relevant_atoms.guess_bonds()  # Guess bonds only for relevant hydrogens


            ag_protein = u.select_atoms(f"(protein {custom_sel}) and ((name H* and bonded (name N* or name O* or name P* or name S*)) or name N* or name O* or name P* or name S*)", updating=True)

        ag_misc = u.select_atoms(f'not (protein or {water})', updating=True) #Keeping this for non-biological systems or where other solvent is important

    except:
        #Make water only
        print('No protein found, creating a network of only waters')
        ag_wat = u.select_atoms(f"resname HOH or resname WAT or resname SOL")

    #extract coordinates from frame of interest
    u.trajectory[frame_idx] 

    #Initiate active site reference atomgroup
    if active_region_reference is not None:
        active_region_residue = u.select_atoms(active_region_reference, updating=True)
    else:
        active_region_residue = None


    #Create network instance
    water_network = WaterNetwork()

    #Check which network type is desired
    if network_type == 'water-protein':
        water_only = False
        #Add protein atoms to network
        for atm in ag_protein.atoms:
            try:
                msa_resid = msa_indexing[atm.resid-1] #CHECK THIS 
            except:
                msa_resid = None
            water_network.add_atom(atm.index, atm.name, atm.resname, *atm.position, atm.resid, msa_resid)
    elif network_type == 'water-water':
        water_only = True
    else:
        raise ValueError("Provide a valid network type. Current valid network types include 'water-protein' or 'water-water'")

    #Add waters to network
    for mol in ag_wat.residues:
        ats = [atom for atom in mol.atoms]
        #Water molecules are objects which contain H1, H2, O atoms
        water_network.add_water(mol.resid, *ats, mol.resid)
    #Either find connections among only oxygens in waters or add hydrogens as well
    if directed:
        water_network.generate_directed_network(u.dimensions, msa_indexing, active_region_residue, active_region_COM=active_region_COM, active_region_radius=active_region_radius, 
                                                active_region_only=active_region_only, water_only=water_only, angle_criteria=angle_criteria, 
                                                max_connection_distance=max_connection_distance)
    else:
        water_network.generate_oxygen_network(u.dimensions, msa_indexing, active_region_residue, active_region_COM=active_region_COM, active_region_radius=active_region_radius, 
                                              active_region_only=active_region_only, water_only=water_only, max_connection_distance=max_connection_distance)
    
    return water_network


def initialize_network(topology_file, trajectory_file, structure_directory='.', network_type='water-protein', 
                       include_hydrogens=False, custom_selection=None, active_region_reference=None, active_region_COM=False, active_region_only=False, 
                       active_region_radius=8.0, water_name=None, multi_model_pdb=False, max_distance=3.0, angle_criteria=None,
                       analysis_conditions='all', analysis_selection='all', project_networks=False, return_network=False, 
                       cluster_coordinates=False, clustering_method='hdbscan', min_cluster_samples=15, eps=None, msa_indexing=True, 
                       alignment_file='alignment.txt', combined_fasta='all_seqs.fa', fasta_directory='fasta', classify_water=False,
                       classification_file_base='DYNAMIC', MSA_reference_pdb=None, water_reference_resids=None,  num_workers=4, shortest_path_nodes=None):

    """
    Initialize and compute all water networks per frame for a trajectory.

    Parameters
    ----------
    structure_directory : str
        Path to the directory containing structure files.
    topology_file : str or None, optional
        Path to the topology file. Default is None.
    trajectory_file : str or None, optional
        Path to the trajectory file. Default is None.
    network_type : {'water-water', 'water-protein'}, optional
        Type of network to construct. Default is 'water-protein'.
    include_hydrogens : bool, optional
        If True, includes hydrogens in the network construction. Default is False.
    custom_selection : str or None, optional
        MDAnalysis selection string for defining a custom protein selection. Default is None.
    active_region_reference : str or None, optional
        MDAnalysis selection string defining the center of the active site. Default is None.
    active_region_COM : bool, optional  
        Whether to take center of mass of active site references, or combine for sphere selection
    active_region_only : bool, optional
        If True, only includes active site atoms in the computed network. This significantly 
        reduces computational cost, especially for directed networks. Default is False.
    active_region_radius : float, optional
        Radius (in Å) defining the active site. Default is 8.0.
    water_name : str or None, optional
        Name of water molecules in the system. Default is None.
    multi_model_pdb : bool, optional
        If True, interprets the input PDB file as a multi-model PDB. Default is False.
    max_distance : float, optional
        Maximum distance (in Å) for defining network connections. Default is 3.0.
    angle_criteria : float or None, optional
        Angle cutoff criteria for hydrogen bonding. Default is None.
    analysis_conditions : {'all', 'active_region', 'not_active_region'}, optional
        Specifies which subset of the system to analyze. Default is 'all'.
    analysis_selection : {'all', 'active_region', 'not_active_region'}, optional
        Specifies the selection criteria for the network analysis. Default is 'all'.
    project_networks : bool, optional
        If True, generates PyMOL visualization files of the computed networks. Default is False.
    return_network : bool, optional
        If True, returns the computed network. Default is False.
    cluster_coordinates : bool, optional
        If True, performs clustering analysis on the network. Default is False.
    clustering_method : {'hdbscan', 'dbscan', 'kmeans'}, optional
        Clustering method to use if clustering is enabled. Default is 'hdbscan'.
    cluster_water_only : bool, optional
        If True, clusters only water molecules, excluding protein atoms. Default is True.
    min_cluster_samples : int, optional
        Minimum number of points per cluster for clustering algorithms. Default is 15.
    eps : float or None, optional
        Epsilon parameter for DBSCAN clustering. Default is None.
    msa_indexing : bool, optional
        If True, performs MSA (multiple sequence alignment) indexing. Default is True.
    alignment_file : str, optional
        Path to the alignment file. If the file does not exist, it will be created. Default is 'alignment.txt'.
    combined_fasta : str, optional
        Name of the combined FASTA file. Default is 'all_seqs.fa'.
    fasta_directory : str, optional
        Directory containing individual FASTA files. Default is 'fasta'.
    classify_water : bool, optional
        If True, classifies water molecules based on MSA indexing. Default is True.
    classification_file_base : str, optional
        Name of outputted csv for classification file. Default is DYNAMIC
    MSA_reference_pdb : str or None, optional
        Path to the reference PDB file for MSA-based selections. Default is None.
    water_reference_resids : list or None, optional
        List of residue indices used for water angle analysis in a specific PDB file. Default is None.
    num_workers : int, optional
        Number of CPU cores to use for parallel computation. Default is 4.
    shortest_path_nodes : list, optional
        List of tuples of nodes to perform shortest path analysis among. Default is None (shortest path among entire network will be returned)

    Returns
    -------
    tuple
        - metrics: dict
            Dictionary of calculated metircs
        - networks: list, optional 
            List of WaterNetwork objects (if return_network=True)
        - centers: dict
            Dictionary of cluster centers (if clustering is on)

    """
    def process_frame(frame_idx, coords=None, ref_coords=None, residues=None, references=None):
        """
        Internal function to make parallelizing each frame easier

        Returns:
        Calculated metrics for given frame
        """

        print(f"Processing frame {frame_idx}")

        #If an MSA has been performed
        if msa_indexing == True:

            #Assuming fasta file is named similarly to the pdb -- need sequence files for MSA alignment
            try:
                fasta_individual = [f for f in os.listdir(fasta_directory) if (topology_file.split('.')[0].split('_')[0] in f and 'fa' in f)][0]
                #Generate MSA if file does not exist and output MSA indices corresponding to partcicular sequence
                msa_indices = sequence_processing.generate_msa_alignment(alignment_file, combined_fasta, os.path.join(fasta_directory, fasta_individual))
            #If MSA cannot be done, use residues as msa_indices
            except:
                print(f'Warning: Could not find an equivalent fasta file for {pdb_file}. Check your naming schemes!')
                msa_indices = residues

        else:
            msa_indices = None

        if active_region_reference is not None and MSA_reference_pdb is not None:
            u = mda.Universe(os.path.join(pdb_dir, pdb_file))
            resids = u.residues.resids.tolist()
            reference_resids, msa_indices_reference = references

            #print(active_region_reference.split())
            #numbers = [f for f in active_region_reference if f.isnumeric()]
            #if len(numbers) == 1:
            #    msa_active_region_ref = sequence_processing.convert_msa_to_individual(msa_indices=msa_indices, msa_indices_ref=msa_indices_reference, resids=resids, resid_sequence_ref=reference_resids, resid_individual_ref=int(active_region_reference.split()[1]))
            #    msa_active_region_ref = f"resid {msa_active_region_ref} {' '.join(active_region_reference.split()[2:])}"
            #else:
            #    msa_numbers = []
            #    for num in numbers:
            #        msa_active_region_num = sequence_processing.convert_msa_to_individual(msa_indices=msa_indices,  msa_indices_ref=msa_indices_reference, resids=resids, resid_sequence_ref=reference_resids, resid_individual_ref=num)
            #        msa_numbers.append(msa_active_region_num)
                    

            # Extract full numbers from the reference string
            numbers = [num for num in active_region_reference.split() if num.isnumeric()]

            if len(numbers) == 1:
                msa_active_region_ref = sequence_processing.convert_msa_to_individual(
                    msa_indices=msa_indices,
                    msa_indices_ref=msa_indices_reference,
                    resids=resids,
                    resid_sequence_ref=reference_resids,
                    resid_individual_ref=int(numbers[0])  # Convert to integer
                )
                msa_active_region_ref = f"resid {msa_active_region_ref} {' '.join(active_region_reference.split()[2:])}"

            else:
                msa_numbers = []
                for num in numbers:
                    msa_active_region_num = sequence_processing.convert_msa_to_individual(
                        msa_indices=msa_indices,
                        msa_indices_ref=msa_indices_reference,
                        resids=resids,
                        resid_sequence_ref=reference_resids,
                        resid_individual_ref=int(num)  # Convert to integer
                    )
                    msa_numbers.append(str(msa_active_region_num))  # Ensure it's a string

                # Replace original numbers with their MSA-mapped equivalents
                words = active_region_reference.split()
                for i, word in enumerate(words):
                    if word in numbers:  # If the word was a number in the original reference
                        words[i] = msa_numbers.pop(0)  # Replace in order

                msa_active_region_ref = " ".join(words)  # Reconstruct the modified string

        else:
            msa_active_region_ref = active_region_reference

        #Create WaterNetwork object
        network = extract_objects_per_frame(pdb_file, traj_file, frame_idx, network_type, custom_selection, active_region_reference, active_region_COM,
                                            active_region_radius=active_region_radius, water_name=water_name, msa_indexing=msa_indices, 
                                            active_region_only=active_region_only, directed=include_hydrogens, angle_criteria=angle_criteria, 
                                            max_connection_distance=max_distance)

        metrics = {}
        #Calculate metrics as per user input
        if analysis_conditions['density'] == 'on':
            metrics['density'] = network.get_density(selection=analysis_selection)
        if analysis_conditions['connected_components'] == 'on':
            metrics['connected_components'] = network.get_connected_components(selection=analysis_selection)
        if analysis_conditions['interaction_counts'] == 'on':
            metrics['interaction_counts'] = network.get_interactions()
        if analysis_conditions['per_residue_interactions'] == 'on':
            metrics['per_residue_interaction'] = network.get_per_residue_interactions(selection=analysis_selection)
        if analysis_conditions['characteristic_path_length'] == 'on':
            metrics['characteristic_path_length'] = network.get_CPL(selection=analysis_selection)
        if analysis_conditions['graph_entropy'] == 'on':
            metrics['entropy'] = network.get_entropy(selection=analysis_selection)
        if analysis_conditions['clustering_coefficient'] == 'on':
            metrics['clustering_coefficient'] = network.get_clustering_coefficient(selection=analysis_selection)

        if analysis_conditions['shortest_path'] == 'on':
            if shortest_path_nodes is None:
                metrics['shortest_path'] = network.get_shortest_path(selection=analysis_selection)
            else:
                metrics['shortest_path'] = []
                for (source, target) in shortest_path_nodes:
                    metrics['shortest_path'].append(network.get_shortest_path(selection=analysis_conditions, source=source, target=target))
        #clustering coefficient -- https://www.annualreviews.org/content/journals/10.1146/annurev-physchem-050317-020915

        #Classify waters
        if classify_water: 
            if msa_indices is None:
                print('ERROR: No MSA indices found, waters cannot be classified without a common indexing reference!')
                print('Did you set msa_indexing=on?')
                raise ValueError
            
            #Select reference coords
            if len(ref_coords) == 1:
                ref2_coords = None
            else:
                ref2_coords = ref_coords[1]

            classification_dict = residue_analysis.classify_waters(network, ref1_coords=ref_coords[0], ref2_coords=ref2_coords)

            #Write classification dict into a csv file
            with open(f'msa_classification/{classification_file_base}.csv', 'a') as FILE:
                for key, val in classification_dict.items():
                        FILE.write(f"{frame_idx},{key},{val[0]},{val[1]}\n")

        #Save coodinates for clustering
        if coords is not None:
            if active_region_only == True:
                selection = 'active_region'
            else:
                selection='all'

            coords = network.get_all_coordinates(selection=selection)
            metrics['coordinates'] = np.array(coords).reshape(-1,3)

        #Create pymol projections for each frame
        if project_networks:
            import WatCon.visualize_structures as visualize_structures
            visualize_structures.pymol_project_oxygen_network(network, filename=f'{frame_idx+2}.pml', out_path='pymol_projections', active_region_only=active_region_only)
        #Do not do this for large trajectories
        if return_network:
            return (metrics, network)
        else:
            return(metrics, None)
        

    
    #Get pdb and traj file
    pdb_dir = structure_directory
    pdb_file = os.path.join(structure_directory, topology_file)
    traj_file = os.path.join(structure_directory, trajectory_file)

    if analysis_conditions == 'all':
        analysis_conditions = {
            'density': 'on',
            'connected_components': 'on',
            'interaction_counts': 'on',
            'per_residue_interactions': 'on',
            'characteristic_path_length': 'on',
            'graph_entropy': 'on',
            'clustering_coefficient': 'on',
            'shortest_path': 'on'
        }

    #Create universe object just once to get number of frames
    try:
        u = mda.Universe(pdb_file, traj_file)
        frames = len(u.trajectory)
        residues = u.residues.resids.tolist()

    except:
        if multi_model_pdb == True:
            u = mda.Universe(pdb_file, multiframe=True)
            frames = len(u.trajectory)
        else:
            print('Warning: You are attempting to create networks for only one structure. Consider using generate_static_networks instead')
            u = mda.Universe(pdb_file)
            frames = 0

    #Initialize empty list to collect coordinates if clustering
    if cluster_coordinates:
        coords = []
    else:
        coords=None

    #Initialize ref_coords if classifying water
    ref_coords = [None]
    if active_region_reference is not None:
        if MSA_reference_pdb is not None:
            u = mda.Universe(os.path.join(pdb_dir, MSA_reference_pdb))
            reference_resids = u.residues.resids.tolist()
            try:
                fasta_individual = [f for f in os.listdir(fasta_directory) if (MSA_reference_pdb.split('_')[0] in f and 'fa' in f)][0] #THIS IS NOT GENERALIZED
            except:
                print(f'Could not find an equivalent fasta file for {MSA_reference_pdb}')

            #Generate MSA alignment if file does not exist and output MSA indices corresponding to partcicular sequence
            msa_indices_reference = sequence_processing.generate_msa_alignment(alignment_file, combined_fasta, os.path.join(fasta_directory, fasta_individual))
            references = [reference_resids, msa_indices_reference]
        else:
            references=None

    if classify_water:
        #Find ref_coords if particular residue is indicated
        if water_reference_resids is not None:
            u = mda.Universe(os.path.join(structure_directory, MSA_reference_pdb))

            #Allows for maximum 2 reference resids
            if isinstance(water_reference_resids, list):
                ref_coords = [u.select_atoms(f"resid {water_reference_resids[0]} and name CA").positions, u.select_atoms(f"resid {water_reference_resids[1]} and name CA").positions]
            else:
                ref_coords = [u.select_atoms(f"resid {water_reference_resids} and name CA").positions]

        #Write header for classification file
        os.makedirs('msa_classification', exist_ok=True)
        with open(f'msa_classification/{classification_file_base}.csv', 'w') as FILE:
            FILE.write('Frame Index,Resid,MSA_Resid,Index_1,Index_2,Protein_Atom,Classification,Protein_Coords,Water_Coords,Angle_1,Angle_2\n')


    #Parallelized so there is one worker allocated for each frame
    results = Parallel(n_jobs=num_workers)(delayed(process_frame)(frame_idx, coords, ref_coords, residues, references) for frame_idx in range(frames))
    network_metrics, networks = zip(*results)
    #Cluster coordinates after networks are created returns metrics and centers
    if cluster_coordinates:
        print('Clustering...')

        # Assuming network_metrics contains dictionaries with 'coordinates' arrays
        coordinates = [
            f['coordinates'] for f in network_metrics if f['coordinates'].shape[1] == 3
        ]
        combined_coordinates = np.concatenate([arr for arr in coordinates], axis=0)

        cluster_labels, cluster_centers = get_clusters(combined_coordinates, cluster=clustering_method, min_samples=min_cluster_samples, eps=eps, n_jobs=num_workers, filename_base=classification_file_base)
        return (network_metrics, networks, cluster_centers)

    #Return only metrics if no clustering
    return(network_metrics, networks, None)
