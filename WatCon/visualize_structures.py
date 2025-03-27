'''
Create pdbs and pml files for visualization of water networks
'''

import os, sys
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

def pymol_project_oxygen_network(network, filename='STATE.pml', out_path='pymol_projections', active_region_only=False, water_only=False):
    """
    Generate a PyMOL (.pml) file to visualize the calculated network.

    Parameters
    ----------
    network : WaterNetwork
        The WaterNetwork object representing the network.
    filename : str
        Name of the output .pml file.
    out_path : str
        Directory where the output file will be saved.
    active_region_only : bool, optional
        Whether to include only active site connections. Default is False.
    water_only : bool, optional
        Whether to show only water-water connections. Default is False.

    Returns
    -------
    None
    """

    os.makedirs(out_path, exist_ok=True)

    with open(os.path.join(out_path, filename), 'w') as FILE:
        #for i, mol in enumerate(network.water_molecules):
            #FILE.write(f"show spheres, id {mol.O.index+1}\nset sphere_scale, 0.3, id {mol.O.index+1}\n")

        #If active_region_only, then only project connections within the active site. Note this is redundant if you have 
        #chosen to only include active site atoms in your whole network
        if active_region_only:
            connection_list = [f for f in network.connections if f[4]=='active_region']
        else:
            connection_list = network.connections

        if water_only:
            connection_list = [f for f in connection_list if f[3]=='WAT-WAT']


        for i, connection in enumerate(connection_list):
            FILE.write(f"distance interaction{i}, id {connection[0]+1}, id {connection[1]+1}\n")
            FILE.write(f"show spheres, id {connection[0]+1}\nshow spheres, id {connection[1]+1}\n")

        #Stylistic preferences
        FILE.write(f"set dash_radius, 0.15, interaction*\nset dash_color, black, interaction*\n")
        FILE.write(f"set dash_gap, 0.0, interaction*\nhide labels, interaction*\n")
        FILE.write(f"set sphere_scale, 0.2\n")
        FILE.write(f"bg white\n")

        #Create WaterNetwork group
        FILE.write(f"group WaterNetwork, interaction*\n")


def project_clusters(coordinate_list, filename_base='CLUSTER',b_factors=None):
    """
    Generate a PDB file to visualize cluster centers.

    Parameters
    ----------
    coordinate_list : array-like, dict
        List of coordinates representing cluster centers.
    filename_base : str
        Naming scheme to use for the outputted PDB file.
    separate_files : bool, optional
        Whether to create separate XYZ files for each cluster or a single combined PDB (default is True).
    b_factors : array-like, optional
        Optional list of values to replace the B-factor column.

    Returns
    -------
    None
    """
    
    #Obtain coordinates from cluster centers
    if type(coordinate_list) == dict:
        coordinate_list = coordinate_list.values() 

    if b_factors is None:
        b_factors = [0.00] * len(coordinate_list)  # Default B-factor is 0.00 if not provided

    os.makedirs('cluster_pdbs', exist_ok=True)

    filename = f"cluster_pdbs/{filename_base}.pdb"
    with open(filename, 'w') as FILE:
        atom_serial = 1
        for label, (center, b_factor) in enumerate(zip(coordinate_list, b_factors)):
            # Write each cluster center as an oxygen atom in PDB format
            FILE.write(
                f"ATOM  {atom_serial:5d}  O   HOH A{label:4d}    "
                f"{center[0]:8.3f}{center[1]:8.3f}{center[2]:8.3f}  1.00 {b_factor:6.2f}           O\n"
            )
            atom_serial += 1


def export_graph_to_pdb(graph, output_file):
    """
    Generate a PDB file with dummy oxygen atoms at the coordinates of a given graph.

    Parameters
    ----------
    graph : networkx.Graph
        The NetworkX graph containing node coordinates.
    output_file : str
        Name of the output PDB file.

    Returns
    -------
    None
    """

    if not output_file.endswith('.pdb'):
        output_file = f"{output_file}.pdb"

    os.makedirs('graph_pdbs', exist_ok=True)
    with open(f"graph_pdbs/{output_file}", 'w') as f:
        atom_serial = 1
        # Collect nodes involved in active site edges
        active_region_nodes = set()
        for edge1, edge2, data in graph.edges(data=True):
            if data.get('active_region') == 'active_region':
                active_region_nodes.add(edge1)
                active_region_nodes.add(edge2)

        # Write nodes as PDB atoms
        for node, data in graph.nodes(data=True):
            if node in active_region_nodes:  # Check if node is part of active site
                x, y, z = data.get('pos', (0.0, 0.0, 0.0))  # Default position if not provided
                f.write(
                    f"ATOM  {atom_serial:5d}  O   HOH A{atom_serial:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           O\n"
                )
                atom_serial += 1

def save_xyz(network, filename='STATE.xyz'):
    """
    Generate an XYZ file containing graph coordinates.

    Parameters
    ----------
    network : WaterNetwork
        The WaterNetwork object containing graph data.
    filename : str
        Name of the output XYZ file.

    Returns
    -------
    None
    """
    if not output_file.endswith('.xyz'):
        output_file = f"{output_file}.xyz"

    with open(filename, 'w') as FILE:
        FILE.write(f"{len(network.molecules)}\n")
        for i, mol in enumerate(network.molecules):
            x, y, z  = mol.coordinates
            FILE.write(f"{mol.name}\t{x:.3f} {y:.3f} {z:.3f}\n")

