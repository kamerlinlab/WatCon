'''
Cluster water coordinates and analyze conservation to clustered networks
'''

import os
import numpy as np
import networkx as nx
from sklearn.cluster import OPTICS, DBSCAN, HDBSCAN
from sklearn.preprocessing import StandardScaler, MinMaxScaler
import pandas as pd
import pickle


def combine_graphs(list_of_graphs):
    """
    Combine multiple NetworkX graph objects into a single graph.

    Parameters
    ----------
    list_of_graphs : list of networkx.Graph
        List of NetworkX graph objects to be merged.

    Returns
    -------
    networkx.Graph
        A combined graph containing all nodes and edges from input graphs.
    """
    U = nx.disjoint_union_all(list_of_graphs)
    return U

def collect_coordinates(pkl_list):
    """
    Collect coordinates into one array from .pkl files

    Parameters
    ----------
    pkl_list : list
        List of pkl_files (full paths)

    Returns
    -------
    np.ndarray
        Array of combined coordinates
    """
    combined_coords = []
    for file in pkl_list:
        with open(file, 'rb') as FILE:
            e = pickle.load(FILE)

        try:
            combined_coords.append([f['coordinates'] for f in e if f['coordinates'].shape[1] == 3])
        except: 
            print('Could not find ["coordinates"], check your inputs!')

    return np.array(combined_coords)

def get_coordinates_from_topology(pdb_file, atom_selection='all'):
    """
    Collect coordinates from a given MDAnalysis-readable topology file

    Parameters
    ----------
    pdb_file : str
        Full path to MDAnalysis-readable topology file

    Returns
    -------
    np.ndarray
        Array of coordinates
    """

    u = mda.Universe(pdb_file)
    ag = u.select_atoms(atom_selection)

    coords = ag.positions
    return np.array(coords).reshape(-1,3)

def cluster_nodes(combined_graph, cluster='hdbscan', min_samples=10):
    """
    Cluster node positions from a combined NetworkX graph.

    Parameters
    ----------
    combined_graph : networkx.Graph
        A combined graph used for clustering.
    cluster : str
        Clustering method, can be 'optics', 'dbscan', or 'hdbscan'.
    min_samples : int
        Minimum number of samples required for a cluster.

    Returns
    -------
    tuple
        - Cluster labels (array-like)
        - Cluster centers (dict)
    """
    node_positions = nx.get_node_attributes(combined_graph, 'pos')
    positions = [pos for (node, pos) in node_positions.items()]
    if cluster == 'optics':
        print('Using OPTICS clustering')
        clustering = OPTICS(max_eps=1.0, metric='euclidean', min_samples=min_samples).fit(positions)
    elif cluster == 'dbscan':
        print('Using DBSCAN clustering')
        clustering = DBSCAN(eps=0.1, min_samples=min_samples).fit(positions)
    elif cluster == 'hdbscan':
        print('Using HDBSCAN Clustering')
        clustering = HDBSCAN(min_cluster_size=min_samples, cluster_selection_epsilon=0.0, algorithm='kd_tree').fit(positions)
    cluster_labels = clustering.labels_
    unique_labels = np.unique(cluster_labels)
    unique_labels.sort() 
    cluster_centers = {}
    for label in unique_labels:
        if label != -1:
            cluster_indices = np.where(cluster_labels == label)[0]
            cluster_center = np.mean(np.array([positions[i] for i in cluster_indices]), axis=0)
    
            cluster_centers[label] = cluster_center
    print(len(cluster_centers))
    return(cluster_labels, cluster_centers)


def cluster_coordinates_only(coordinate_list, cluster='hdbscan', min_samples=10, eps=0.0, n_jobs=1):
    """
    Cluster a set of coordinates.

    Parameters
    ----------
    coordinate_list : list of array-like
        A combined list of all coordinates to be clustered.
    cluster : str
        Clustering method, can be 'optics', 'dbscan', or 'hdbscan'.
    min_samples : int
        Minimum number of samples required for a cluster.

    Returns
    -------
    tuple
        - Cluster labels (array-like)
        - Cluster centers (dict)
    """
    try:
        coordinate_list = np.array(coordinate_list).reshape(-1,3)
    except:
        print("Couldn't reshape coordinates correctly, check your inputs.")

    #scaler = MinMaxScaler()
    #scaler.fit(coordinate_list)
    #coordinate_norm = scaler.transform(coordinate_list)

    if cluster == 'optics':
        print('Using OPTICS clustering')
        clustering = OPTICS(min_samples=min_samples, eps=eps, n_jobs=n_jobs).fit(coordinate_list)
    elif cluster == 'dbscan':
        print('Using DBSCAN clustering')
        clustering = DBSCAN(min_samples=min_samples, eps=eps, n_jobs=n_jobs).fit(coordinate_list)
    elif cluster == 'hdbscan':
        print('Using HDBSCAN clustering')
        print(min_samples, eps, coordinate_list.shape)
        clustering = HDBSCAN(min_cluster_size=min_samples, cluster_selection_epsilon=eps, algorithm='kd_tree', n_jobs=n_jobs).fit(coordinate_list)

    cluster_labels = clustering.labels_
    unique_labels = np.unique(cluster_labels)
    unique_labels.sort() 
    cluster_centers = {}
    for label in unique_labels:
        if label != -1:
            cluster_indices = np.where((cluster_labels == label) & (cluster_labels!=-1))[0]
            #cluster_center = [coordinate_list[i] for i in cluster_indices][0]
            cluster_center = np.mean(np.array([coordinate_list[i] for i in cluster_indices]), axis=0)
    
            cluster_centers[label] = cluster_center

    print(len(cluster_centers))
    return(cluster_labels, cluster_centers)


def find_commonality(networks, centers, names):
    """
    Find the commonality of a list of networks relative to a summary network created from clustering.

    Parameters
    ----------
    networks : list of WaterNetwork
        List of WaterNetwork objects to be analyzed.
    centers : array-like
        Locations of clustered centers.

    Returns
    -------
    dict
        A dictionary containing calculated commonalities for each network.
    """

    commonality_dict = {}
    dist = lambda x1, y1, z1, x2, y2, z2: np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
    for i, network in enumerate(networks):
        conserved = 0
        unique = 0
        net = networks[i]
        for wat in net.water_molecules:
            x1 = wat.O.coordinates[0]
            y1 = wat.O.coordinates[1]
            z1 = wat.O.coordinates[2]
            if any((dist(x1,y1,z1, x2,y2,z2)<3) for (x2, y2, z2) in centers):
                local_waters_count = len([wat for wat in net.water_molecules if (dist(wat.O.coordinates[0], wat.O.coordinates[1],wat.O.coordinates[2],x1,y1,z1)<6 and dist(wat.O.coordinates[0], wat.O.coordinates[1],wat.O.coordinates[2],x1,y1,z1)>2)])+1
                conserved += 1/local_waters_count
            else:
                unique += 1
        commonality_dict[names[i]] = conserved/len(centers)
    return commonality_dict


def identify_conserved_water_clusters(networks, centers, dist_cutoff=1.0, filename_base='CLUSTERS'):
    """
    Create a dictionary of cluster conservation and generate a PDB file of clusters.

    Parameters
    ----------
    networks : list of WaterNetwork
        List of WaterNetwork objects to analyze.
    centers : array-like
        List of XYZ coordinates representing cluster centers.
    dist_cutoff : float
        Distance cutoff to classify a water molecule as part of a cluster.
    filename_base : str
        Base filename for saving projected clusters.

    Returns
    -------
    dict
        A dictionary mapping cluster centers to the count of included waters.
    """
    from sklearn.preprocessing import MinMaxScaler
    from WatCon.visualize_structures import project_clusters
    #Add protein atoms eventually
    center_dict = {}
    dist = lambda x1, y1, z1, x2, y2, z2: np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
    for i, network in enumerate(networks):
        for j, center in enumerate(centers):
            if str(j) not in center_dict.keys():
                center_dict[str(j)] = 0
            for wat in network.water_molecules:
                x1 = wat.O.coordinates[0]
                y1 = wat.O.coordinates[1]
                z1 = wat.O.coordinates[2]
                x2, y2, z2 = center
                if dist(x1,y1,z1, x2,y2,z2) <= dist_cutoff:
                    center_dict[str(j)] += 1

    values = np.array([f for f in center_dict.values()]).reshape(-1,1)
    scaler = MinMaxScaler()
    scaler.fit(values)

    b_factors = scaler.transform(values).flatten()
    project_clusters(centers, filename_base=filename_base, b_factors=b_factors, separate_files=False)

    return(center_dict)


def create_clustered_network(clusters, max_connection_distance, create_graph=True):
    """
    Create a WaterNetwork object from cluster centers.

    Parameters
    ----------
    clusters : array-like
        List of XYZ coordinates of cluster centers.
    max_connection_distance : float
        Maximum allowed distance between two clusters to form an interaction.
    create_graph : bool, optional
        Whether to create a NetworkX graph from the clustered WaterNetwork. Default is True.

    Returns
    -------
    WaterNetwork
        A WaterNetwork object representing the clustered water network.
    """
    from WatCon.generate_static_networks import WaterNetwork, WaterAtom, WaterMolecule

    clustered_network = WaterNetwork()

    for i, center in enumerate(clusters):
        o = WaterAtom(i, 'O', i, *center)

        water = WaterMolecule(i, o, H1=None, H2=None, residue_number=i)
        clustered_network.water_molecules.append(water)

    clustered_network.connections = clustered_network.find_connections(dist_cutoff=max_connection_distance, water_only=True)

    if create_graph:
        G = nx.Graph()
        for molecule in clustered_network.water_molecules:
                G.add_node(molecule.O.index, pos=molecule.O.coordinates, atom_category='WAT', MSA=None) #have nodes on all oxygens
        for connection in clustered_network.connections:
                G.add_edge(connection[0], connection[1], connection_type=connection[3], active_region=connection[4])

        clustered_network.graph = G
    return clustered_network

def identify_conserved_water_interactions_clustering(networks, clusters, max_connection_distance=2.0, dist_cutoff=1.0, filename_base='CLUSTER'):
    """
    Rank water-water interactions in relation to clustering.

    Parameters
    ----------
    networks : list of WaterNetwork
        List of WaterNetwork objects to be analyzed.
    centers : array-like
        List of XYZ coordinates of cluster centers.
    max_connection_distance : float
        Maximum allowed distance between two clusters to form an interaction.
    dist_cutoff : float
        Distance cutoff to classify a water molecule as part of a cluster.
    filename_base : str
        Base filename for saving projected clusters.

    Returns
    -------
    dict
        A dictionary mapping cluster interactions to the count of included waters.
    """
    from sklearn.preprocessing import MinMaxScaler
    import matplotlib.pyplot as plt

    interaction_dict = {}
    dist = lambda x1, y1, z1, x2, y2, z2: np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)

    clustered_network = create_clustered_network(clusters, max_connection_distance, create_graph=True)

    for network in networks:
        for connection in clustered_network.connections:
            name = f"{connection[0]}, {connection[1]}"
            if str(name) not in interaction_dict.keys():
                interaction_dict[str(name)] = 0
            
            x2, y2, z2 = clusters[connection[0]]
            x3, y3, z3 = clusters[connection[1]]

            if any(dist(wat.O.coordinates[0], wat.O.coordinates[1], wat.O.coordinates[2], x2,y2,z2) <= dist_cutoff for wat in network.water_molecules) and any(dist(wat.O.coordinates[0], wat.O.coordinates[1], wat.O.coordinates[2], x3,y3,z3) <= dist_cutoff for wat in network.water_molecules):
                interaction_dict[str(name)] += 1

    values = np.array([f for f in interaction_dict.values()]).reshape(-1,1)
    scaler = MinMaxScaler()
    scaler.fit(values)

    b_factors = scaler.transform(values).flatten()    

    pairs = [f for f in interaction_dict.keys()]
    num_interactions = len(pairs)

    # Sort pairs based on b_factors
    pairs = [pair for _, pair in sorted(zip(b_factors, pairs), key=lambda x: x[0])]

    cmap = plt.get_cmap('bwr')

    len_colors = np.linspace(0,1,num_interactions)
    colors = cmap(len_colors)

    #colors = [(int(r*255), int(g*255), int(b*255)) for r, g, b, _ in colors]
    colors = [(float(r), float(g), float(b)) for r, g, b, _ in colors]
    colors = [color for _, color in sorted(zip(b_factors, colors), key=lambda x: x[0])]


    with open(f'{filename_base}.pml', 'w') as f:
        for i, (pair, color) in enumerate(zip(pairs, colors)):
            f.write(f"distance interaction{i}, resid {pair.split(',')[0]}, resid {pair.split(',')[1]}\n")
            f.write(f"set dash_color, [{color[0]},{color[1]},{color[2]}], interaction{i}\n")
            f.write(f"show spheres, resid {pair.split(',')[0]}\n")
            f.write(f"show spheres, resid {pair.split(',')[1]}\n")
        
        f.write("hide labels, all\n")
        f.write('set dash_radius, 0.15, interaction*\n')    
        f.write('set dash_gap, 0.0, interaction*\n')
        f.write('hide labels, interaction*\n')
        f.write('set sphere_scale, 0.2\n')
        f.write('bg white\n')
        f.write('group ClusterNetwork, interaction*\n')

    return interaction_dict


def identify_clustered_angles(classification_file, ref1_coords, ref2_coords):
    """
    NOTE: NOT TESTED YET

    MAKE IT SO THAT YOU CLUSTER FOR CRYSTAL STRUCTURES AND TAKE MINIMUM COORDINATES FOR MD
    """
    from scipy.optimize import minimize

    '''
    def angle_constraint(wat,prot,ref,theta):
        wat = np.array(wat)
        prot = np.array(prot)
        ref = np.array(ref)

        v1 = wat - prot
        v2 = ref - prot

        cos_theta = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        return np.cos(theta) - cos_theta
    '''

    def angle_constraint(wat, *args):
        prot, ref, theta = args
        return np.cos(theta) - np.dot(wat - prot, ref - prot) / (np.linalg.norm(wat - prot) * np.linalg.norm(ref - prot))


    def dist_norm(wat, prot):
        return np.linalg.norm(wat - prot)
    
    def find_wat_coords(prot_coords, ref1_coords, ref2_coords, theta1, theta2):
        theta1 = np.radians(theta1)
        theta2 = np.radians(theta2)

        w0 = (np.array(ref1_coords)+np.array(ref2_coords))/2
        w0 = (np.array(ref1_coords) + np.array(ref2_coords)) / 2 + 0.1 * (np.array(prot_coords) - w0)

        constraints = (
            {'type': 'eq', 'fun': angle_constraint, 'args': (prot_coords, ref1_coords, theta1)},
            {'type': 'eq', 'fun': angle_constraint, 'args': (prot_coords, ref2_coords, theta2)}
        )

        result = minimize(dist_norm, w0, args=(prot_coords,), constraints=constraints, method='trust-constr')

        if result.success:
            return result.x
        else:
            raise ValueError('Optimization failed')

    df = pd.read_csv(classification_file, delimiter=',')

    classification_dict = {}
    coord_dict = {}

    for i, row in df.iterrows():
        if str(row['MSA_Resid']) not in classification_dict.keys():
            classification_dict[str(row['MSA_Resid'])] = []
            coord_dict[str(row['MSA_Resid'])] = []
        classification_dict[str(row['MSA_Resid'])].append((float(row['Angle_1']), float(row['Angle_2'])))
        coord_dict[str(row['MSA_Resid'])].append(np.array([float(f) for f in row['Protein_Coords'].split()]))


    cluster_conservation_dict = {}

    for msa_resid, values in classification_dict.items():
        if str(msa_resid) not in cluster_conservation_dict.keys():
            cluster_conservation_dict[str(msa_resid)] = {}

        hdb = HDBSCAN(min_cluster_size=3)
        values = np.array([np.array(f) for f in values]).reshape(-1,2)

        if values.shape[0] > 3:
            hdb.fit(values)

            cluster_labels = hdb.labels_
            unique_labels = np.unique(cluster_labels)
            unique_labels.sort() 
            cluster_centers = {}


            for label in [f for f in unique_labels if f != -1]:
                cluster_indices = np.where(cluster_labels == label)[0]
                cluster_center = np.mean(np.array([values[i] for i in cluster_indices]), axis=0)
                cluster_centers[label] = cluster_center

            centroid_coords = {}
            for label, center in cluster_centers.items():
                min_distance=1000
                for val_ind, val in enumerate(values):
                    if np.linalg.norm([center, val]) < min_distance:
                        closest_coord = coord_dict[str(msa_resid)][val_ind]

                centroid_coords[label] = closest_coord

            for label in cluster_labels:
                if label != -1:
                    if label not in cluster_conservation_dict[str(msa_resid)].keys():
                        cluster_conservation_dict[str(msa_resid)][str(label)] = {'counts': 0, 'center':cluster_centers[label], 'closest_coord': centroid_coords[label]}

                        print(centroid_coords[label], ref1_coords, ref2_coords, *cluster_centers[label])
                        wat_coords = find_wat_coords(centroid_coords[label], ref1_coords, ref2_coords, *cluster_centers[label])
                        cluster_conservation_dict[str(msa_resid)][str(label)] = {'wat_coord': wat_coords}
                    cluster_conservation_dict[str(msa_resid)][str(label)]['counts'] += 1

    return cluster_conservation_dict


def plot_consevation_angles(cluster_conservation_dict, output_filebase='angle_clusters', output_dir='pymol_projections'):
    with open(os.path.join(output_dir, f"{output_filebase}.pml")) as FILE:
        for msa in cluster_conservation_dict.items():
            for i, label in enumerate(cluster_conservation_dict[msa].keys()):
                water_coord = cluster_conservation_dict[msa][label]['wat_coord']
                prot_coord = cluster_conservation_dict[msa][label]['closest_coord']

                strength = cluster_conservation_dict[msa][label]['wat_coord']


                FILE.write(f"pseudoatom {msa}_{i}_protein, pos={prot_coord}\n")
                FILE.write(f"pseudoatom {msa}_{i}_water, pos={water_coord}\n")
                FILE.write(f"distance interaction_{msa}_{i}, {msa}_{i}_protein, {msa}_{i}_water\n")
        
        FILE.write("hide labels, all\n")
        FILE.write('set dash_radius, 0.15, interaction*\n')    
        FILE.write('set dash_gap, 0.0, interaction*\n')
        FILE.write('hide labels, interaction*\n')
        FILE.write('set sphere_scale, 0.2\n')
        FILE.write('set sphere_color, oxygen, *_water')
        FILE.write('bg white\n')
        FILE.write('group AngleInteractions, interaction*\n')
        FILE.write('group PseudoProteins, *_protein\n')
        FILE.write('group PseudoWater, *_water\n')




def find_clusters_from_densities(density_file, output_name=None, threshold=1.5):
    """
    Find clusters from densities.

    Parameters
    ----------
    density_file : str
        .dx file containing density information
    output_name : str
        Base name for output
    threshold : float
        Threshold for cutoff

    Returns
    -------
    np.ndarray
        Array of density hotspot coordinate locations


    """
    from gridData import Grid
    import scipy.ndimage as ndimage
    if output_name is None:
        output_name = f"{density_file.split('.dx')[0]}"

    grid = Grid(density_file)
    data = grid.grid
    origin = np.array(grid.origin)
    delta = np.array(grid.delta)

    neighborhood = np.ones((3,3,3))
    local_max = (data == ndimage.maximum_filter(data, footprint=neighborhood))
    hotspot_indices = np.argwhere(local_max & (data > threshold))

    hotspot_coords = np.array([origin + idx * delta] for idx in hotspot_indices)

    with open(f"{output_name}.pdb", 'w') as FILE:
        for i, (x,y,z) in enumerate(hotspot_coords, start=1):
            FILE.write(f"ATOM{i:5d}  O   HOH     1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           O  \n")

    return(hotspot_coords)

    
def plot_commonality(files, input_directory, cluster_pdb, plot_type='bar', output='commonality'):
    """
    Plot commonality to cluster centers.

    Parameters
    ----------
    files : list
        List of outputted .pkl files from WatCon
    input_directory : str
        Directory containing files
    cluster_pdb : str
        PDB of clusters
    plot_type : {'bar', 'hist'}
        Type of plot
    output : str
        Base filename for saving images

    Returns
    -------
    None
    """
    import matplotlib.pyplot as plt

    name_list = []
    network_list = []

    names = True
    #Check if files have names (from crystal structures) or do not
    with open(os.path.join(input_directory,files[0]), 'rb') as FILE:
        e = pickle.load(FILE)
        if len(e) < 4:
            names = False
        
    if names:
        for i, file in enumerate(files):
            with open(os.path.join(input_directory,file), 'rb') as FILE:
                e = pickle.load(FILE)
            name_list.extend(e[3])
            network_list.extend(e[1])
    else:
        for i, file in enumerate(files):
            with open(os.path.join(input_directory,file), 'rb') as FILE:
                e = pickle.load(FILE)

            name_list.extend([f"{i}-{j}" for j, _ in enumerate(e)])
            network_list.extend(e[1])

    with open(cluster_pdb, 'r') as FILE:
        lines = FILE.readlines()

    centers = np.zeros((len(lines), 3))
    for i, line in enumerate(lines):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            centers[i,:] = x,y,z
    
    commonality_dict = find_commonality(network_list, centers, name_list)

    if plot_type == 'bar':
        fig, ax = plt.subplots(1, figsize=(5,3), tight_layout=True)
        gene_data = {}
        for i, (name, commonality) in enumerate(commonality_dict.items()):

            # Plot the bar graph
            names = list(commonality_dict.keys())
            x = np.arange(len(names))
            width = 0.5

        for i, name in enumerate(names):
            ax.bar(x[i], commonality, width, color='gray', hatch='//', edgecolor='k', label=name)

        ax.set_xticks(x)
        ax.set_xticklabels(names, fontsize=12)
        ax.set_ylabel('Conservation score', fontsize=15)
        ax.tick_params(axis='y', labelsize=12)
        ax.tick_params(axis='x', rotation=60)
        #ax.legend(frameon=True, edgecolor='k', fontsize=10) 
        
        plt.savefig(f"{output}_bar.png", dpi=200)
    
    elif plot_type == 'hist':
        fig, ax = plt.subplots(1, figsize=(3,2), tight_layout=True)
        vals = np.array(list(commonality_dict.values()))

        hist, xedges = np.histogram(vals, density=True, bins=15)
        xcenters = (xedges[1:]+xedges[:-1])/2

        ax.plot(xcenters, hist)
        ax.set_xlabel('Commonality score')
        ax.set_ylabel('Density')
        plt.savefig(f"{output}_hist.png", dpi=200)
        
    else:
        print('Select a valid plot type. Currently only "bar" or "hist".')
        raise ValueError