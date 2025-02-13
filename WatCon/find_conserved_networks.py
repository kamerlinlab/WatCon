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


def cluster_nodes(combined_graph, cluster='hdbscan', min_samples=10, eps=0.6):
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
        clustering = OPTICS(min_samples=min_samples).fit(positions)
    elif cluster == 'dbscan':
        print('Using DBSCAN clustering')
        clustering = DBSCAN(eps=eps, min_samples=min_samples).fit(positions)
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


def cluster_coordinates_only(coordinate_list, cluster='hdbscan', min_samples=10, n_jobs=1, eps=0.6):
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
        clustering = OPTICS(min_samples=min_samples, n_jobs=n_jobs).fit(coordinate_list)
    elif cluster == 'dbscan':
        print('Using DBSCAN clustering')
        clustering = DBSCAN(eps=eps, min_samples=min_samples, n_jobs=n_jobs).fit(coordinate_list)
    elif cluster == 'hdbscan':
        print('Using HDBSCAN clustering')
        clustering = HDBSCAN(min_cluster_size=min_samples, cluster_selection_epsilon=0.0, algorithm='kd_tree', n_jobs=n_jobs).fit(coordinate_list)

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


def project_clusters(cluster_centers, filename_base='CLUSTER', separate_files=True):
    """
    Create pdb file to visualize cluster centers

    Parameters:
    - cluster_centers: Cluster centers outputted from clustering function
    - filename_base: Naming scheme to use for outputted pdb
    - separate_files (True/False): Indicate whether to create a separate xyz for each 
                                   cluster or one combined pdb

    Returns:
    None
    """
    if separate_files == True:
        for label, center in cluster_centers.items():
            #Create a filename for the cluster
            filename = f"{filename_base}_{label}.xyz"
            with open(filename, 'w') as FILE:
                # Write the cluster center's coordinates to the file
                FILE.write(f"O\t{center[0]:.3f} {center[1]:.3f} {center[2]:.3f}\n")
    else:
        filename = f"{filename_base}.pdb"
        with open(filename, 'w') as FILE:
            atom_serial = 1
            for label, center in cluster_centers.items():
                # Write each cluster center as an oxygen atom in PDB format
                FILE.write(
                    f"ATOM  {atom_serial:5d}  O   HOH A{label:4d}    "
                    f"{center[0]:8.3f}{center[1]:8.3f}{center[2]:8.3f}  1.00  0.00           O\n"
                )
                atom_serial += 1

def find_commonality(networks, centers):
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
    for network in networks.keys():
        conserved = 0
        unique = 0
        net = networks[network]
        for wat in net.water_molecules:
            x1 = wat.coordinates[0]
            y1 = wat.coordinates[1]
            z1 = wat.coordinates[2]
            if any((dist(x1,y1,z1, x2,y2,z2)<3) for (x2, y2, z2) in centers.values()):
                local_waters_count = len([wat for wat in net.water_molecules if (dist(wat.coordinates[0], wat.coordinates[1],wat.coordinates[2],x1,y1,z1)<6 and dist(wat.coordinates[0], wat.coordinates[1],wat.coordinates[2],x1,y1,z1)>2)])+1
                conserved += 1/local_waters_count
            else:
                unique += 1
        commonality_dict[network] = conserved/len(centers.keys())
    return commonality_dict


if __name__ == "__main__":
    pdbs = os.listdir('aligned_pdbs_with_water')
    pdbs.sort()
    find_conserved_networks_static('aligned_pdbs_with_water', pdbs)
    #find_conserved_networks_MD('/home/abrownless3/Documents/Summer2024/get_PTP1B_traject_short/analysis/output_fixed.gro', '/home/abrownless3/Documents/Summer2024/get_PTP1B_traject_short/analysis/output_final.xtc', file_format='gro')
    

'''
#OLD 

def find_conserved_networks_static(pdb_dir, pdbs, fasta_dir='fasta',alignment_file=None,combined_fasta=None, water_only=True, optional_id=None):
    networks = {}
    graphs = []
    densities = []
    for pdb in pdbs:
        if not water_only:
            try:
                fasta = [f for f in os.listdir(fasta_dir) if pdb.split('_aligned')[0] in f][0]
                networks[pdb.split('.')[0]] = initialize_network(os.path.join(pdb_dir, pdb), alignment_file, combined_fasta, os.path.join(fasta_dir, fasta)) 
            except FileNotFoundError as e:
                print('Fasta file not found.')
        else:
             networks[pdb.split('.')[0]] = initialize_network(os.path.join(pdb_dir, pdb)) 
        G = networks[pdb.split('.')[0]].generate_network()
        density = networks[pdb.split('.')[0]].density()
        densities.append(density)
        #pymol_project_network_old(networks[pdb.split('.')[0]], f"{pdb.split('.')[0]}.pml")
        graphs.append(G)

    combined_graph = combine_graphs(graphs)
    labels, centers = cluster_nodes(combined_graph, cluster='dbscan', min_samples=3)
    #project_clusters(labels, combined_graph, optional_id)
    #project_summary_network(centers, optional_id)
    commonality_dict = find_commonality(networks, centers)
    return(commonality_dict, densities)

def find_conserved_networks_MD(positions, cluster='optics', min_samples=10):
    if cluster == 'optics':
        clustering = OPTICS(min_samples=4).fit(positions)
    elif cluster == 'dbscan':
        clustering = DBSCAN(eps=1.0, min_samples=min_samples).fit(positions)

    cluster_labels = clustering.labels_
    unique_labels = np.unique(cluster_labels)
    unique_labels.sort() 
    cluster_centers = {}
    for label in unique_labels[1:]:
        cluster_indices = np.where(cluster_labels == label)[0]
        cluster_center = np.mean(np.array([positions[i] for i in cluster_indices]), axis=0)

        cluster_centers[label] = cluster_center
    return(cluster_labels, cluster_centers)


def create_predicted_network(cluster_centers, positions, out_name='SUMMARY'):
    dist = lambda x1, y1, z1, x2, y2, z2: np.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
    with open(f"{out_name}.pdb", 'w') as FILE:
        atom_index = 1
        residue_index = 1
        for coord1 in cluster_centers:
            for coord2 in positions:
                if dist(*coord1, *coord2) > 1:
                            # Write nodes as PDB atoms
                    FILE.write(
                        f"ATOM  {atom_index:5d}  O   HOH A{residue_index:4d}    {coord1[0]:8.3f}{coord1[1]:8.3f}{coord1[2]:8.3f}  1.00  0.00           O\n"
                    )
                    atom_index += 1
                    residue_index += 1
'''
