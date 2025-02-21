'''
Cluster water coordinates and analyze conservation to clustered networks
'''

import os
import numpy as np
import networkx as nx
#from  generate_dynamic_networks import initialize_water_network, pymol_project_oxygen_network, pymol_project_directed_network
from generate_static_networks import initialize_network
from sequence_processing import generate_msa_alignment
from sklearn.cluster import OPTICS, DBSCAN, HDBSCAN
from sklearn.preprocessing import StandardScaler, MinMaxScaler


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
    for label in unique_labels[1:]:
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

