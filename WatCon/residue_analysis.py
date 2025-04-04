'''
Per-residue water analysis
'''

import os
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import distances
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd


def get_interaction_counts(network, selection='all'):
    """
    Calculate numbers of interactions split by 'water-water' and 'water-protein'.

    Parameters
    ----------
    network : WaterNetwork object
    selection : {'all', 'active_region', 'not_active_region'}
        Specifies which subset of the graph to analyze.

    Returns
    ----------
    dict
        Describes number of 'water-water' and 'water-protein' interactions.
    """
    interaction_counts = {'water-water': 0, 'water-protein': 0}
    for _, _, data in network.graph.edges(data=True):
        if selection=='all':
            if data['connection_type']=='WAT-PROT':
                interaction_counts['water-protein'] += 1
            else:
                interaction_counts['water-water'] += 1
        else:
            if data['active_region'] == selection:
                if data['connection_type']=='WAT-PROT':
                    interaction_counts['water-protein'] += 1
                else:
                    interaction_counts['water-water'] += 1
    return interaction_counts


def get_per_residue_interactions(network, selection='all', msa=False):
    """
    Calculate numbers of interactions per residue.

    Parameters
    ----------
    network : WaterNetwork object
    selection : {'all', 'active_region', 'not_active_region'}
        Specifies which subset of the graph to analyze.
    msa : bool, optional
        Indicate whether to use msa common residue indices

    Returns
    ----------
    dict
        Describes number of interactions per residue.
    """
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
        if data['connection_type'] == 'WAT-PROT' and (selection == 'all' or data['active_region'] == selection)
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
    Classify all water-protein interactions based on two reference angles.

    This function analyzes the geometric relationships between water molecules and 
    protein atoms by calculating angles relative to two reference points.

    Parameters
    ----------
    network : WaterNetwork
        The water network object containing interaction data.
    ref1_coords : array-like
        Coordinates of the first reference point (e.g., an atom or centroid).
    ref2_coords : array-like
        Coordinates of the second reference point (e.g., an atom or centroid).

    Returns
    -------
    dict
        A dictionary describing interaction classifications and calculated angles.
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

    #print('ref1:', ref1_coords)
    #print('ref2:', ref2_coords)
    
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
        
    if np.all(ref1_coords == ref2_coords): 
        print('Reference coordinates are the same, input two separate coordinates.')
        raise ValueError
        #return None
    
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
            prot_name = [f"{f.resid},{f.msa_resid},{connection[0]},{connection[1]},{connection[2]},{connection[5]},{f.coordinates[0]} {f.coordinates[1]} {f.coordinates[2]},{wat_coords[0]} {wat_coords[1]} {wat_coords[2]}" for f in network.protein_atoms if (f.index == connection[0] or f.index == connection[1])][0]
            classification_dict[prot_name] = [angle1, angle2] #Consider combining into one value
        else:
            continue

    return classification_dict
    

def plot_interactions_from_angles(csvs, input_dir='msa_classification',output_dir='MSA_images'):
    """
    Plot water classifications from 2-angle analysis

    Parameters
    ----------
    csvs: list
        List of .csv files outputted from classify_waters
    output_dir: str
        Directory to write images

    Returns
    ---------
    None

    """

    dfs = {}
    for csv in csvs:
        df = pd.read_csv(os.path.join(input_dir,csv), delimiter=',')
        name = '_'.join(csv.split('_')[0:2]).split('.')[0]
        dfs[name] = df

    scatters = {}
    classifications = {}
    max_length=0
    MSA_min = 1000
    MSA_max = 1

    pdb_names = {}
    for name, df in dfs.items():
        if len([f for f in df.iterrows()]) > max_length:
            max_length = len([f for f in df.iterrows()])
        df.sort_values(by='MSA_Resid')
        scatters[name] = {}
        classifications[name] = {}
        pdb_names[name] = {}
        for i, row in df.iterrows():
            if row['MSA_Resid'] in scatters[name].keys():
                classifications[name][row['MSA_Resid']].append(row['Classification'])
                scatters[name][row['MSA_Resid']].append((row['Angle_1'], row['Angle_2']))
                pdb_names[name][row['MSA_Resid']].append(row['PDB ID'])
            else:
                scatters[name][row['MSA_Resid']] = [(row['Angle_1'], row['Angle_2'])]
                classifications[name][row['MSA_Resid']] = [row['Classification']]
                pdb_names[name][row['MSA_Resid']] = [row['PDB ID']]
            if row['MSA_Resid'] < MSA_min:
                MSA_min = row['MSA_Resid']
            if row['MSA_Resid'] > MSA_max:
                MSA_max = row['MSA_Resid']

    # Collect all unique MSAs across all scatters
    all_MSAs = sorted({MSA for data in scatters.values() for MSA in data.keys()})

    # Calculate global x and y limits
    all_x = []
    all_y = []
    for data in scatters.values():
        for coord_list in data.values():
            for coords in coord_list:
                x, y = map(float, coords)
                all_x.append(x)
                all_y.append(y)

    x_min, x_max = min(all_x), max(all_x)
    y_min, y_max = min(all_y), max(all_y)

    # Add padding for better visualization
    x_range = x_max - x_min
    y_range = y_max - y_min
    x_min, x_max = x_min - 0.1 * x_range, x_max + 0.1 * x_range
    y_min, y_max = y_min - 0.1 * y_range, y_max + 0.1 * y_range



    names = list(scatters.keys())
    print(names)

    colors = {'DYNAMIC': 'gray', 'STATIC': {'backbone': 'dodgerblue', 'sidechain': 'mediumorchid'}}  # Adjust colors
    os.makedirs(output_dir, exist_ok=True)

    # Generate and save a separate plot for each MSA
    for MSA in all_MSAs:
        plt.figure(figsize=(3, 2.5), tight_layout=True)

        # Collect all MD data
        md_x, md_y = [], []
        
        for name, data in scatters.items():
            if MSA in data:
                if name.startswith('DYNAMIC'):  # Combine all MD_* data
                    x_vals, y_vals = zip(*data[MSA])  # Extract coordinates
                    md_x.extend(map(float, x_vals))
                    md_y.extend(map(float, y_vals))

        # Plot the combined MD data as a surface
        if md_x and md_y:
            hist, x_edges, y_edges = np.histogram2d(md_x, md_y, bins=50)
            X, Y = np.meshgrid(x_edges[:-1], y_edges[:-1])
            mesh = plt.contourf(X, Y, np.log(hist.T), cmap="gray")
            cbar = plt.colorbar(mesh)
            cbar.set_label('log(Density)')

        # Plot STATIC data as scatter plots
        for name, data in scatters.items():
            if MSA in data and name == 'STATIC':
                for i, coords in enumerate(data[MSA]):
                    x, y = map(float, coords)
                    classification = classifications[name][MSA][i]
                    color = colors['STATIC']['backbone'] if 'backbone' in classification else colors['STATIC']['sidechain']
                    name_new = pdb_names[name][MSA][i]
                    if 'open' in name_new or 'Open' in name_new:
                        facecolor='none'
                    else:
                        facecolor=color
                    plt.scatter(x, y, edgecolor=color, facecolor=facecolor, s=10)
                    
                    #plt.text(x,y, name_new)

        plt.xlim(x_min, x_max)
        plt.ylim(y_min, y_max)
        plt.xticks(np.arange(0,181, 50))
        plt.yticks(np.arange(0,181,50))
        plt.xlabel('Angle 1')
        plt.ylabel('Angle 2')
        plt.title(f"Common residue {int(MSA)}", fontsize=11)
        #plt.title(f"MSA: {int(MSA)}", fontsize=12)
        plt.tight_layout()
        #plt.savefig(f"{output_dir}/MSA_{int(MSA)}.png", dpi=200)
        plt.savefig(f"{output_dir}/MSA_{int(MSA)}.png", dpi=600)
        plt.close()


def histogram_metrics(all_files, input_directory, concatenate):
    """
    Plot histograms for calculated metrics

    Parameters
    ----------
    all_files: list
        List of all files
    input_directory: str
        Directory which contains .pkl files
    concatenate: list
        List of files to concatenate

    Returns
    ---------
    None

    """
    import pickle

    #Initialize dictionaries to store data
    metrics = ['density', 'characteristic_path_length', 'entropy']

    metric_dict = {'density':[],'characteristic_path_length':[], 'entropy':[], 'water-water':[], 'water-protein':[]}

    metrics_plot =  ['density', 'characteristic_path_length', 'entropy', 'water-water', 'water-protein']

    #Formatted titles for plotting
    plotting_names = ['Graph Density', 'CPL', 'Graph Entropy', 'Water-Water', 'Water-Protein']


    #Combine data in all concatenated files
    for file in concatenate:
        watcon_file = os.path.join(input_directory, file)
        with open(watcon_file, 'rb') as FILE:
            e = pickle.load(FILE)

        for ts_dict in e[0]:
            metric_dict['water-water'].append(ts_dict['interaction_counts']['water-water'])
            metric_dict['water-protein'].append(ts_dict['interaction_counts']['water-protein'])

            for metric in metrics:
                if isinstance(ts_dict[metric], float):
                    metric_dict[metric].append(ts_dict[metric])
                else:
                    metric_dict[metric].extend([f for f in ts_dict[metric]])

    #Select all other files
    all_files = [f for f in all_files if f not in concatenate]

    #Create list to store other dictionaries
    metric_dicts = []
    for file in all_files:
        metric_dict_static = {'density':[],'characteristic_path_length':[], 'entropy':[], 'water-water':[], 'water-protein':[]}
        watcon_file = os.path.join(input_directory, file)
        with open(watcon_file, 'rb') as FILE:
            e = pickle.load(FILE)

        for ts_dict in e[0]:
            metric_dict_static['water-water'].append(ts_dict['interaction_counts']['water-water'])
            metric_dict_static['water-protein'].append(ts_dict['interaction_counts']['water-protein'])
            for metric in metrics:
                if isinstance(ts_dict[metric], float):
                    metric_dict_static[metric].append(ts_dict[metric])
                else:
                    metric_dict_static[metric].extend(f for f in ts_dict[metric])
        metric_dicts.append(metric_dict_static)

    #Begin plotting
    for i, metric in enumerate(metrics_plot):
        fig, ax = plt.subplots(1,figsize=(3,2))
        fig.subplots_adjust(left=0, right=0.85)
        metric_cur_concatenate = np.array(metric_dict[metric])


        hist, xedges = np.histogram(metric_cur_concatenate, bins=15, density=True)
        xcenters = (xedges[1:]+xedges[:-1])/2
        ax.plot(xcenters, hist, label='Concatenated')

        #sns.kdeplot(data=np.array(metric_cur_dynamic), ax=ax, bw_adjust=2)


        for i, metric_dict in enumerate(metric_dicts):
            metric_cur = np.array(metric_dict[metric])


            hist, xedges = np.histogram(metric_cur, bins=15, density=True)
            xcenters = (xedges[1:]+xedges[:-1])/2
            ax.plot(xcenters, hist, label=f'Sample {i}')

            #sns.kdeplot(data=np.array(metric_cur_static), ax=ax)

        ax.legend(fontsize=8, frameon=False)
        ax.set_xlabel(plotting_names[i])
        ax.set_ylabel('Density')

        fig.savefig(f"{metric}_comparehists.png", dpi=200, bbox_inches='tight')