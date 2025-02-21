'''
Create pdbs and pml files for visualization of water networks
'''

import os, sys
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

def pymol_project_oxygen_network(network, filename='STATE.pml', out_path='pymol_projections', active_site_only=False, water_only=False):
    """
    Create .pml file to showcase calculated network

    Parameters:
    - network: WaterNetwork object
    - filename: Output filename
    - out_path: Directory for output
    - active_site_only: Indicate whether to include only active site connections or all connections
    - water_only: Indiacte whether to show connections among waters only

    Returns:
    None
    """
    with open(os.path.join(out_path, filename), 'w') as FILE:
        #for i, mol in enumerate(network.water_molecules):
            #FILE.write(f"show spheres, id {mol.O.index+1}\nset sphere_scale, 0.3, id {mol.O.index+1}\n")

        #If active_site_only, then only project connections within the active site. Note this is redundant if you have 
        #chosen to only include active site atoms in your whole network
        if active_site_only:
            connection_list = [f for f in network.connections if f[4]=='active_site']
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


def project_clusters(coordinate_list, filename_base='CLUSTER', separate_files=True):
    """
    Create pdb file to visualize cluster centers

    Parameters:
    - cluster_centers: List of coordinates or cluster centers
    - filename_base: Naming scheme to use for outputted pdb
    - separate_files (True/False): Indicate whether to create a separate xyz for each 
                                   cluster or one combined pdb

    Returns:
    None
    """
    
    #Obtain coordinates from cluster centers
    if type(coordinate_list) == dict:
        coordinate_list = coordinate_list.values() 

    if separate_files == True:
        for label, center in enumerate(coordinate_list):
            #Create a filename for the cluster
            filename = f"{filename_base}_{label}.xyz"
            with open(filename, 'w') as FILE:
                # Write the cluster center's coordinates to the file
                FILE.write(f"O\t{center[0]:.3f} {center[1]:.3f} {center[2]:.3f}\n")
    else:
        filename = f"{filename_base}.pdb"
        with open(filename, 'w') as FILE:
            atom_serial = 1
            for label, center in enumerate(coordinate_list):
                # Write each cluster center as an oxygen atom in PDB format
                FILE.write(
                    f"ATOM  {atom_serial:5d}  O   HOH A{label:4d}    "
                    f"{center[0]:8.3f}{center[1]:8.3f}{center[2]:8.3f}  1.00  0.00           O\n"
                )
                atom_serial += 1


def export_graph_to_pdb(graph, output_file):
    """
    Create a pdb file with dummy oxygen atoms at coordinates within a given graph

    Parameters:
    - graph: networkx graph object
    - output_file: Name of output file

    Returns:  
    None
    """

    if not output_file.endswith('.pdb'):
        output_file = f"{output_file}.pdb"

    with open(output_file, 'w') as f:
        atom_serial = 1
        # Collect nodes involved in active site edges
        active_site_nodes = set()
        for edge1, edge2, data in graph.edges(data=True):
            if data.get('active_site') == 'active_site':
                active_site_nodes.add(edge1)
                active_site_nodes.add(edge2)

        # Write nodes as PDB atoms
        for node, data in graph.nodes(data=True):
            if node in active_site_nodes:  # Check if node is part of active site
                x, y, z = data.get('pos', (0.0, 0.0, 0.0))  # Default position if not provided
                f.write(
                    f"ATOM  {atom_serial:5d}  O   HOH A{atom_serial:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           O\n"
                )
                atom_serial += 1

def save_xyz(network, filename='STATE.xyz'):
    """
    Create xyz file with graph coordinates

    Parameters:
    - network: WaterNetwork object
    - filename: Name of output file
    """
    if not output_file.endswith('.xyz'):
        output_file = f"{output_file}.xyz"

    with open(filename, 'w') as FILE:
        FILE.write(f"{len(network.molecules)}\n")
        for i, mol in enumerate(network.molecules):
            x, y, z  = mol.coordinates
            FILE.write(f"{mol.name}\t{x:.3f} {y:.3f} {z:.3f}\n")

