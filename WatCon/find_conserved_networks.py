'''
Cluster water coordinates and analyze conservation to clustered networks
'''

import os
import numpy as np
import networkx as nx
from sklearn.cluster import OPTICS, DBSCAN, HDBSCAN
from sklearn.preprocessing import StandardScaler, MinMaxScaler
import pandas as pd


def combine_graphs(list_of_graphs):
    """
    Combine networkx graph objects into one graph

    Parameters:
    - list_of_graphs: List of networkx graph objects

    Returns:
    - Combined graph
    """
    U = nx.disjoint_union_all(list_of_graphs)
    return U


def cluster_nodes(combined_graph, cluster='hdbscan', min_samples=10):
    """
    Cluster nodes positions from combined networkx graph

    Parameters: 
    - combined_graph: Combined graph for clustering
    - cluster: Clustering method, can be 'optics', 'dbscan', or 'hdbscan'
    - min_samples: Minimum samples for cluster

    Returns:
    - Cluster labels, cluster centers
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
    Cluster coordinates only

    Parameters:
    - coordinate_list: Combined list of all coordinates
    - cluster: Clustering method, can be 'optics', 'dbscan', or 'hdbscan'
    - min_samples: Minimum samples for cluster

    Returns:
    - Cluster labels, cluster centers
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
    Find commonality of a list of networks to a summary network created from clustering

    Parameters:
    - networks: List of WaterNetwork objects
    - centers: Locations of clustered centers

    Returns:
    - Dictionary of calculated commonalities for each network 
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
    Create dictionary of conservation of clusters and create pdb of clusters

    Parameters:
    - networks: List of WaterNetwork object
    - centers: List of xyz coordinates of cluster centers
    - dist_cutoff: Distance cutoff to classify a water as within a cluster
    - filename_base: Name to save the projected clusters

    Returns:
    - Dictionary of cluster centers with counts of included waters
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
    project_clusters(centers, filename_base=filename_base, b_factors=b_factors)

    return(center_dict)


def create_clustered_network(clusters, max_connection_distance, create_graph=True):
    """
    Create WaterNetwork object from cluster centers

    Parameters:
    - clusters: List of xyz coordinates of clusters
    - max_connection_distance: Maximum distance between two clusters to interact
    - create_graph: Create a networkx object from the clustered WaterNetwork

    Returns;
    - WaterNetwork object
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
                G.add_edge(connection[0], connection[1], connection_type=connection[3], active_site=connection[4])

        clustered_network.graph = G
    return clustered_network

def identify_conserved_water_interactions_clustering(networks, clusters, max_connection_distance=3.0, dist_cutoff=1.0, filename_base='CLUSTER'):
    """
    Create dictionary which ranks water-water interactions in relation to clustering

    Parameters:
    - networks: List of WaterNetwork objects
    - centers: List of xyz coordinates of cluster centers
    - max_connection_distance: Maximum distance between two clusters to interact
    - dist_cutoff: Distance cutoff to classify a water as within a cluster
    - filename_base: Name to save the projected clusters

    Returns:
    - Dictionary of cluster interactions with counts of included waters
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


def identify_conserved_waterprotein_interactions_angles(classification_file):
    df = pd.read_csv(classification_file, delimiter=',')

    classification_dict = {}

    for i, row in df.iterrows():
        if str(row['MSA_Resid']) not in classification_dict.keys():
            classification_dict[str(row['MSA_Resid'])] = []
        classification_dict[str(row['MSA_Resid'])].append((float(row['Angle_1']), float(row['Angle_2'])))


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

            for label in cluster_labels:
                if label != -1:
                    if label not in cluster_conservation_dict[str(msa_resid)].keys():
                        cluster_conservation_dict[str(msa_resid)][str(label)] = {'counts': 0, 'center':cluster_centers[label]}
                    cluster_conservation_dict[str(msa_resid)][str(label)]['counts'] += 1

    return cluster_conservation_dict