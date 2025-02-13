'''
Per-residue water analysis
'''

import os
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import distances
import networkx as nx
import matplotlib.pyplot as plt

#from generate_dynamic_networks import *

'''
def get_waterprotein_interactions(network_group):
    # Initialize output dictionaries and lists
    interaction_data = {
        'No-active-site': {
            'Water-Protein': [[], []],  # [waterprotein_int_list, all_watprot]
            'Water-Water': [[], []]     # [watwat_int_list, all_watwat]
        },
        'Active-site': {
            'Water-Protein': [[], []],  # [activeprotein_list, all_activeprot]
            'Water-Water': [[], []]     # [activewater_list, all_activewat]
        }
    }

    # Process each network in the group
    for network in network_group:
        # Counters and connection lists
        waterprotein_interactions = 0
        water_interactions = 0
        activeprotein = 0
        activewater = 0
        
        waterprotein_connections = []
        waterwater_connections = []
        activeprotein_connections = []
        activewater_connections = []

        # Iterate through edges in the network graph
        for cn1, cn2, data in network.graph.edges(data=True):
            if data['connection_type'] == 'WAT-PROT':
                waterprotein_interactions += 1
                waterprotein_connections.append((cn1, cn2))

                if data['active_site'] == 'active_site':
                    activeprotein += 1
                    activeprotein_connections.append((cn1, cn2))
            else:
                water_interactions += 1
                waterwater_connections.append((cn1, cn2))

                if data['active_site'] == 'active_site':
                    activewater += 1
                    activewater_connections.append((cn1, cn2))

        # Append interactions to the relevant lists
        interaction_data['No-active-site']['Water-Protein'][0].append(waterprotein_interactions)
        interaction_data['No-active-site']['Water-Protein'][1].append(waterprotein_connections)

        interaction_data['No-active-site']['Water-Water'][0].append(water_interactions)
        interaction_data['No-active-site']['Water-Water'][1].append(waterwater_connections)

        interaction_data['Active-site']['Water-Protein'][0].append(activeprotein)
        interaction_data['Active-site']['Water-Protein'][1].append(activeprotein_connections)

        interaction_data['Active-site']['Water-Water'][0].append(activewater)
        interaction_data['Active-site']['Water-Water'][1].append(activewater_connections)

    return interaction_data


def get_waterprotein_interactions_singleframe(network):
    # Initialize output dictionaries and lists
    interaction_data = {
        'No-active-site': {
            'Water-Protein': [],  # [waterprotein_int_list, all_watprot]
            'Water-Water': []    # [watwat_int_list, all_watwat]
        },
        'Active-site': {
            'Water-Protein': [],  # [activeprotein_list, all_activeprot]
            'Water-Water': []    # [activewater_list, all_activewat]
        }
    }

    waterprotein_interactions = 0
    water_interactions = 0
    activeprotein = 0
    activewater = 0
        
    waterprotein_connections = []
    waterwater_connections = []
    activeprotein_connections = []
    activewater_connections = []
    
    # Iterate through edges in the network graph
    for cn1, cn2, data in network.graph.edges(data=True):
        if data['connection_type'] == 'WAT-PROT':
            waterprotein_interactions += 1
            waterprotein_connections.append((cn1, cn2))

            if data['active_site'] == 'active_site':
                activeprotein += 1
                activeprotein_connections.append((cn1, cn2))
        else:
            water_interactions += 1
            waterwater_connections.append((cn1, cn2))

            if data['active_site'] == 'active_site':
                activewater += 1
                activewater_connections.append((cn1, cn2))

    # Append interactions to the relevant list
    interaction_data['No-active-site']['Water-Protein'].append(waterprotein_interactions)
    interaction_data['No-active-site']['Water-Protein'].append(waterprotein_connections)

    interaction_data['No-active-site']['Water-Water'].append(water_interactions)
    interaction_data['No-active-site']['Water-Water'].append(waterwater_connections)

    interaction_data['Active-site']['Water-Protein'].append(activeprotein)
    interaction_data['Active-site']['Water-Protein'].append(activeprotein_connections)

    interaction_data['Active-site']['Water-Water'].append(activewater)
    interaction_data['Active-site']['Water-Water'].append(activewater_connections)

    return interaction_data
'''

def get_interaction_counts(network, selection='all'):
    interaction_counts = {'water-water': 0, 'water-protein': 0}
    for _, _, data in network.graph.edges(data=True):
        if selection=='all':
            if data['connection_type']=='WAT-PROT':
                interaction_counts['water-protein'] += 1
            else:
                interaction_counts['water-water'] += 1
        else:
            if data['active_site'] == selection:
                if data['connection_type']=='WAT-PROT':
                    interaction_counts['water-protein'] += 1
                else:
                    interaction_counts['water-water'] += 1
    return interaction_counts

'''
def get_per_residue_interactions(network, selection='all'):
    residue_dict = {}
    if selection == 'all':
        edges = [(_,_,data) for (_,_,data) in network.graph.edges(data=True) if data['connection_type']=='WAT-PROT']
    else:
        edges = [(_,_,data) for (_,_,data) in network.graph.edges(data=True) if data['active_site']==selection and data['connection_type']=='WAT-PROT']
    for (cn1, cn2, data) in edges:
        try:
            cn1_res = [f.resid for f in network.water_molecules if cn1==f.O.index][0]
            cn1_wat = True
        except:
            cn1_res = [f.resid for f in network.protein_atoms if cn2==f.index][0]
            cn1_wat = False
        if cn1_wat == True:
            cn2_res = [f.resid for f in network.protein_atoms if cn2==f.index][0]
        else:
            cn2_res = [f.resid for f in network.water_molecules if cn2==f.O.index][0]

        if cn1_wat == True:
            if cn2_res in residue_dict.keys():
                residue_dict[str(cn2_res)] += 1
            else:
                residue_dict[str(cn2_res)] = 1
        else:
            if cn1_res in residue_dict.keys():
                residue_dict[str(cn1_res)] += 1
            else:
                residue_dict[str(cn1_res)] = 1
    return residue_dict
'''

def get_per_residue_interactions(network, selection='all', msa=False):
    def get_resid_by_index(index, is_water):
        if is_water:
            matches = [f.resid for f in network.water_molecules if f.O.index == index]
        else:
            if msa==True:
                matches = [f.msa_resid for f in network.protein_atoms if f.index == index]
            else:            
                matches = [f.resid for f in network.protein_atoms if f.index == index]
        return matches[0] if matches else None

    residue_dict = {}
    edges = [
        (cn1, cn2, data) for (cn1, cn2, data) in network.graph.edges(data=True)
        if data['connection_type'] == 'WAT-PROT' and (selection == 'all' or data['active_site'] == selection)
    ]

    for cn1, cn2, _ in edges:
        cn1_res = get_resid_by_index(cn1, is_water=True) or get_resid_by_index(cn1, is_water=False)
        cn2_res = get_resid_by_index(cn2, is_water=False) or get_resid_by_index(cn2, is_water=True)
        cn1_wat = cn1_res in [f.resid for f in network.water_molecules]

        if cn1_wat:
            target_res = cn2_res
        else:
            target_res = cn1_res

        if target_res:
            residue_key = str(target_res)
            residue_dict[residue_key] = residue_dict.get(residue_key, 0) + 1

    return residue_dict

'''
def get_per_residue_interactions(network_group, selection='No-active-site'):
    interaction_data = get_waterprotein_interactions(network_group)
    active_interactions = interaction_data[selection]['Water-Protein']
    residue_dict_main = {}
    residue_dict_main['num_interactions'] = {}
    for mol in network_group[0].protein_atoms:
        residue_dict_main['num_interactions'][str(mol.resid)] = []
    for ts, (n_interactions, connections) in enumerate(zip(active_interactions[0], active_interactions[1])):
        network = network_group[ts]
        residue_dict_tmp = {}
        #CHECK THAT YOU"RE NOT DOUBLE COUNTING
        for connection in connections:
            residue = [mol.resid for mol in network.protein_atoms if ((mol.index == connection[0]) or (mol.index == connection[1]))][0] 
            if str(residue) in residue_dict_tmp.keys():
                residue_dict_tmp[str(residue)] += 1
            else:
                residue_dict_tmp[str(residue)] = 1
        for main_key in residue_dict_main['num_interactions'].keys():
            if main_key in residue_dict_tmp.keys():
                residue_dict_main['num_interactions'][main_key].append(residue_dict_tmp[main_key])
            else:
                residue_dict_main['num_interactions'][main_key].append(0)
    return residue_dict_main, interaction_data
'''

'''
def get_interaction_percentage(network_group):
    interaction_dict = {}
    for network in network_group:
        for (u, v) in network.edges():
            #CHECK THAT YOU'RE NOT DOUBLE COUNTING HERE
            if str(u) not in interaction_dict.keys():
                interaction_dict[str(u)] = 1
            else:
                interaction_dict[str(u)] += 1
            if str(v) not in interaction_dict.keys():
                interaction_dict[str(v)] = 1
            else:
                interaction_dict[str(v)] += 1 
    frames = len(network_group)
    for key, values in interaction_dict.items():
        interaction_dict[key] = values/frames
    return interaction_dict
'''


def get_all_water_distances(network_group, box, selection='No-active-site', msa=False, offset=0):
    tmp_list = []
    residue_dict, interaction_data = get_per_residue_interactions(network_group, selection=selection, msa=msa)
    residue_dict['water_distances'] = {}
    for network, connections in zip(network_group, interaction_data[selection]['Water-Protein'][1]):
        dist_arr = []
        for connection in connections:
            residue_atom = [mol for mol in network.protein_atoms if (mol.index == connection[0]) or (mol.index == connection[1])][0]
            if msa==True:
                tmp_list.append(residue_atom.resid+offset) #Add optional offset 
                res = residue_atom.resid
            else:
                tmp_list.append(residue_atom.msa_resid+offset)
                res = residue_atom.msa_resid
            #ONLY WORKS FOR OXYGEN NETWORK RN
            try:
                water_mol = [mol for mol in network.water_molecules if (mol.O.index == connection[0]) or (mol.O.index == connection[1]) 
                          or (mol.H1.index == connection[0]) or (mol.H1.index == connection[1]) or 
                          (mol.H2.index == connection[0]) or (mol.H2.index == connection[1])][0]
            except:
                print(connection[0], connection[1])
            dist = np.min(distances.distance_array(np.array(residue_atom.coordinates), 
                                                   np.array([np.array(water_mol.O.coordinates), np.array(water_mol.H1.coordinates), np.array(water_mol.H2.coordinates)]).reshape(-1,3), box=box))
            dist_arr.append(dist)
        residue_dict['water_distances'][str(res)] = dist_arr
    print(set(tmp_list))
    return residue_dict, interaction_data


def classify_waters(network, ref1_coords, ref2_coords):
    """
    Classify all water-protein interactions based on two reference angles 

    Parameters:
    - network: WaterNetwork object
    - ref1_coords: Coordinates of first reference point
    - ref2_coords: Coordinates of second reference point

    Returns:
    Dictionary describing interactions and calculated angles
    """
    #Maybe extend to implement ML model to optimize angles

    try:
        if ref1_coords[0] is None:
            ref1_coords = [(0,10,0)]
    except:
        if ref1_coords is None:
            ref1_coords = [(0,10,0)]
            
    if ref2_coords is None:
        ref2_coords = [(10,0,10)]
    
    #THERE IS SOME WEIRD STUFF GOING ON WITH THESE SHAPES
    if type(ref1_coords) == tuple:
        ref1_coords=tuple(ref1_coords[0][0]) #Careful
    else:
        ref1_coords=tuple(ref1_coords[0]) #Careful

    #THERE IS SOME WEIRD STUFF GOING ON WITH THESE SHAPES
    if type(ref2_coords) == tuple:
        ref2_coords=tuple(ref2_coords[0][0]) #Careful
    else:
        ref2_coords=tuple(ref2_coords[0]) #Careful
        
    if np.all(ref1_coords == ref2_coords): #Change this to be a try/except with error returned
        print('Reference coordinates are the same, input two separate coordinates.')
        return None
    
    def get_angles(wat_coords, prot_coords, ref_coords):
        v1 = np.array([prot_coords[0]-wat_coords[0], prot_coords[1]-wat_coords[1], prot_coords[2]-wat_coords[2]])
        v2 = np.array([ref_coords[0]-wat_coords[0], ref_coords[1]-wat_coords[1], ref_coords[2]-wat_coords[2]])

        mag_v1 = np.sqrt((v1[0]**2+v1[1]**2+v1[2]**2))
        mag_v2 = np.sqrt((v2[0]**2+v2[1]**2+v2[2]**2))
        angle = (180/np.pi) * np.arccos(np.dot(v1, v2)/(mag_v1*mag_v2))
        return angle

    classification_dict = {}
    for connection in [f for f in network.connections if f[3]=='WAT-PROT']:
        wat_coords = [f.O.coordinates for f in network.water_molecules if (
            f.O.index == connection[0] or
            f.O.index == connection[1] or
            (f.H1 is not None and (f.H1.index == connection[0] or f.H1.index == connection[1])) or
            (f.H2 is not None and (f.H2.index == connection[0] or f.H2.index == connection[1]))
            )][0]
        
        #Maybe change to using CA coordinates if this is too sensitive
        prot_coords = [f.coordinates for f in network.protein_atoms if (f.index == connection[0] or f.index==connection[1])][0]
        if len(prot_coords) > 0:
            angle1 = get_angles(wat_coords, prot_coords, ref_coords=ref1_coords)
            angle2 = get_angles(wat_coords, prot_coords, ref_coords=ref2_coords)
            prot_name = [f"{f.resid}, {f.msa_resid}, {connection[0]}, {connection[1]}, {connection[2]}, {connection[5]}" for f in network.protein_atoms if (f.index == connection[0] or f.index == connection[1])][0]
            classification_dict[prot_name] = [angle1, angle2] #Consider combining into one value
        else:
            continue

    return classification_dict

def plot_waterprotein_interactions(network_group):
    all_interaction_dict = get_waterprotein_interactions(network_group)
    active_site_interaction = get_waterprotein_interactions(network_group)
    fig, ax = plt.subplots(2, figsize=(6,4), tight_layout=True)
    n_frames = len(network_group)
    frames = np.arange(n_frames, 1)
    ax[0].plot(frames, all_interaction_dict['Water_Protein'], color='k', label='Full protein')
    ax[1].plot(frames, active_site_interaction['Water_Protein'], color='r', label='Active site only')
    fig.supxlabel('Frame')
    fig.supylabel('Number of protein/water interactions')
    plt.show()

if __name__ == '__main__':
    parm = ''
    traj = ''
    
    network_group, box, md_universe = initialize_trajectory_waterprotein_network(parm, traj, custom_selection='resname CSP or resname CYM')




