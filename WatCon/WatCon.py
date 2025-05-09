import os, sys
import argparse
import pickle

def parse_inputs(filename):
    """
    Parse inputs from input file

    Parameters
    ----------
    filename : str
        Path of input file

    Returns
    ----------
    tuple
        - str
            structure type (static/dynamic)
        - dict
            Dictionary of kwargs
    """

    #Initialize blank kwargs dict
    kwargs = {}

    analysis_conditions = {}

    with open(filename, 'r') as FILE:
        lines = FILE.readlines()

    for i, line in enumerate(lines):
        if 'Property calculation' in line:
            analysis_conditions = {}
            for analysis_condition in lines[i+1:]:
                if len(analysis_condition) < 1 or 'analysis_selection' in analysis_condition:
                    break
                else:
                    kw = analysis_condition.split(':')[0]
                    kw_value = analysis_condition.split(':')[1].split()[0]
                    analysis_conditions[kw] = kw_value

        #Skip past comments    
        if line.startswith(';') or ':' not in line:
            continue

        #Choose dynamic or static
        elif 'structure_type' in line:
            structure_type = line.split(':')[1].split()[0]

        elif 'make_fastas' in line:
            make_fastas = line.split(':')[1].split()[0]

        #Assign all other kwargs
        else:
            kw = line.split(':')[0]
            kw_value = ' '.join(line.split(':')[1].split()[0].split("#"))
            if kw_value == 'on':
                kw_value = True
            elif kw_value == 'off':
                kw_value = False

            if kw not in analysis_conditions.keys():
                kwargs[kw] = kw_value

    if len(analysis_conditions.keys()) == 0:
        analysis_conditions = 'all'
    kwargs['analysis_conditions'] = analysis_conditions

    if 'water_reference_resids' in kwargs.keys():
        final_resids = []
        resids = kwargs['water_reference_resids'].split(',')
        for resid in resids:
            final_resids.append(int(resid))
        kwargs['water_reference_resids'] = final_resids

    #Change strings to floats/ints
    if 'max_distance' in kwargs.keys():
        kwargs['max_distance'] = float(kwargs['max_distance'])
    
    if 'angle_criteria' in kwargs.keys():
        if kwargs['angle_criteria'] == 'None':
            kwargs['angle_criteria'] = None
        else:
            kwargs['angle_criteria'] = float(kwargs['angle_criteria'])

    if 'active_region_radius' in kwargs.keys():
        kwargs['active_region_radius'] = float(kwargs['active_region_radius'])

    if 'multi_model_pdb' in kwargs.keys():
        if kwargs['multi_model_pdb'] == 'False':
            kwargs['multi_model_pdb'] = False
        else:
            kwargs['multi_model_pdb'] = True
     
    if 'min_cluster_samples' in kwargs.keys():
        kwargs['min_cluster_samples'] = int(kwargs['min_cluster_samples'])

    if 'eps' in kwargs.keys():
        kwargs['eps'] = float(kwargs['eps'])

    if 'num_workers' in kwargs.keys():
        kwargs['num_workers'] = int(kwargs['num_workers'])
    
    if 'water_name' in kwargs.keys():
        if kwargs['water_name'] == 'default':
            kwargs['water_name'] = None
    
    if 'trajectory_name' in kwargs.keys() and kwargs['trajectory_name'] == 'None':
        kwargs['trajectory_name'] = None

    if 'topology_name' in kwargs.keys() and kwargs['topology_name'] == 'None':
        kwargs['topology_name'] = None

    if 'MSA_reference' in kwargs.keys() and kwargs['MSA_reference'] == 'None':
        kwargs['MSA_reference'] = None
    
    return structure_type, kwargs

#NOTE NEED TO MAKE IT SO THAT YOU CAN TELL IT TO MAKE THE FASTAS FROM THE PDBS

def check_conditions(kwargs):
    print('Checking conditions from input file')
    #Consider making some text to check conditions
    pass
    
def parse_analysis(filename):
    """
    Parse inputs from input file

    Parameters
    ----------
    filename : str
        Path of input file

    Returns
    ----------
    dict
        Dictionary of kwargs
    """
    kwargs = {}

    with open(filename, 'r') as FILE:
        lines = FILE.readlines()

    for line in lines:
        if 'concatemate:' in line:
            files_to_concatenate = list(line.split(":")[1].split(';')[0])
            kwargs['concatenate'] = files_to_concatenate
        elif 'active_region_definition' in line:
            active_region_definition = ' '.join(line.split(":")[1].split(';')[0])
            kwargs['active_region_definition'] = active_region_definition

        else:
            if ':' in line:
                kw = line.split(':')[0]
                kw_value =line.split(':')[1].split()[0]
                if kw_value == 'on':
                    kw_value = True
                elif kw_value == 'off':
                    kw_value = False

                kwargs[kw] = kw_value
    return (kwargs)

def run_watcon(structure_type, kwargs):
    """
    Run WatCon from input file

    Parameters
    ----------
    structure_type : {'static', 'dynamic'}
        Type of input structures
    kwargs : dict
        Key word arguments from parsed input file

    Returns
    ----------
    tuple
        Results from initialize network
    """
    if structure_type == 'static':
        from WatCon.generate_static_networks import initialize_network
    else:
        from WatCon.generate_dynamic_networks import initialize_network

    results = initialize_network(**kwargs)
    return results

def run_watcon_postanalysis(concatenate=None, input_directory='watcon_output', histogram_metrics=False, residue_interactions=False, 
                        reference_topology=None, interaction_cutoff=0.0, calculate_densities=False, density_pdb=None, 
                        traj_directory=None, active_region_definition=None, image_output_dir='images',
                        custom_selection=None, water_name=None, cluster_concatenated=False, cluster_method='hdbscan', eps=0.0,
                        n_jobs=1, min_samples=100, cluster_filebase='CLUSTER', calculate_commonality=None, color_by_conservation=None, 
                        classify_waters=False, csv_dir='msa_classification'):
    """
    Run WatCon analysis from input file

    Parameters
    ----------
    concatenate : list
        List of files to concatenate when performing analysis. Default is None
    input_directory : str
        Directory which contains input files. Default is 'watcon_output'
    histogram_metrics : bool
        Indicate whether to histogram calculated metrics. Default is False.
    residue_interactiosn: bool
        Indicate whether to plot bar graphs of water-protein interactions at the residue level. Default is False.
    reference_topology : str
        Full path to reference topology to be used for residue interaction analysis
    calculate_densities : bool
        Calculate densities for combined simulation data. Default is False.
    density_pdb : str
        PDB to use for density calculations. Default is None.
    traj_directory : str
        Directory containing trajectories to be used when calculating densities. Default is None.
    active_region_definition : str
        MDAnalysis selection language to define active site when calculating densities. Default is None.
    image_output_dir : str
        Output directory for images. Default is 'images'
    custom_selection : str
        MDAnalysis selection language for classifying custom residues as 'protein'. Default is None.
    water_name : str
        Residue name of water (not necessary if water is WAT, SOL, H2O, or HOH). Default is None.
    cluster_concatenated : bool
        Indicate whether to cluster coordinates across concatenate list. Default is False.
    cluster_method : {'hdbscan', 'optics', 'dbscan'}
        Cluster method (if cluster_concatenated). Default is 'hdbscan'
    eps : float
        Eps value for clustering method (if cluster_concatenated). Default is 0.0.
    min_samples : int
        Min sample value for clustering method (if cluster_concatenated). Default is 100.
    n_jobs : int
        Number of available cores for clustering 
    cluster_filebase : str
        Base of output files. Default is 'CLUSTER'. OR Desired clustered pdb to calculate commonality
    calculate_commonality : {'bar', 'hist', None}
        Indicate whether to calculate commonality score to a set of clusters and 
        plot either as a bar graph (helpful for discrete static structures) or histogram (helpful for high dimensional dynamic data)
    color_by_conservation : {'all', 'centers', 'connections', None}
        Color cluster centers by conservation, connections among cluster centers, or neither
    classify_waters : bool
        Indicate whether to plot distributions of 2-angle calculations. Default is False.
    csv_dir: str
        Directory containing csvs (for classify_waters)

    Returns
    ----------
    None

    Notes
    ----------
    calculate_densities is NOT recommended for static structures due to sparsity in water positions

    cluster_concatenated is NOT recommended for collections of large trajectories, calculate_densities is recommended instead

    If f"{cluster_filebase}.pdb" already exists and cluster_concatenated == False, then the current clustering pdb will be used for analysis
    """

    import WatCon.residue_analysis as residue_analysis

    #Find all .pkl files in input_directory
    all_files = [f for f in os.listdir(input_directory) if f.endswith('.pkl')]

    #Histogram metircs
    if histogram_metrics:
        residue_analysis.histogram_metrics(all_files, input_directory, concatenate, image_output_dir)

    if residue_interactions:
        if reference_topology is not None:
            residue_analysis.plot_residue_interactions(reference_topology, interaction_cutoff, watcon_directory=input_directory, output_dir=image_output_dir)
        else:
            print('Please provide a reference topology for residue interaction analysis')

    #Calculate densities and save centers as pdb atoms
    if calculate_densities:
        from WatCon.generate_dynamic_networks import collect_densities

        #Check that pdb and trajectory files have been supplied
        if any([density_pdb,traj_directory,active_region_definition] == None):
            print('If calculating densities, must provide valid PDB and trajectory information')
            raise ValueError
        
        trajectories = [os.listdir(traj_directory)]
        trajectories = [os.path.join(traj_directory, f) for f in trajectories]
        density_pdb = os.path.join(traj_directory, density_pdb)

        #Allow for user-defined water name
        if water_name is None:
            water_name = 'resname HOH or resname WAT or resname SOL or resname H2O'
        else:
            water_name = f"resname {water_name}"

        collect_densities(density_pdb, trajectories, active_region_definition, custom_selection, water_name, f"{cluster_filebase}.dx")

    #Cluster concatenated coordinates
    if cluster_concatenated:
        from WatCon.find_conserved_networks import collect_coordinates, cluster_coordinates_only
        from WatCon.visualize_structures import project_clusters
        files = [os.path.join(input_directory, f) for f in concatenate]
        combined_coordinates = collect_coordinates(files)

        cluster_labels, cluster_centers = cluster_coordinates_only(combined_coordinates, cluster='hdbscan', min_samples=min_samples, eps=eps, n_jobs=n_jobs)
        project_clusters(cluster_centers, filename_base=cluster_filebase, separate_files=False, b_factors=None)
    
    if calculate_commonality:
        from WatCon.find_conserved_networks import plot_commonality

        if not os.path.exists(f"cluster_pdbs/{cluster_filebase}.pdb"):
            print(f"Clustered PDB does not exist. Cannot calculate commonality")
            raise ValueError

        plot_commonality(all_files, input_directory, f"cluster_pdbs{cluster_filebase}.pdb" , commonality_dict=None, plot_type=calculate_commonality)

    if color_by_conservation is not None:
        from WatCon.find_conserved_networks import identify_conserved_water_clusters, identify_conserved_water_interactions_clustering, get_coordinates_from_topology
        
        centers = get_coordinates_from_topology(f"cluster_pdbs/{cluster_filebase}.pdb")

        for file in os.listdir(input_directory):
            watcon_file = os.path.join(input_directory, file)
            with open(watcon_file, 'rb') as FILE:
                e = pickle.load(FILE)
            network_metrics, networks, centers, names = e
            name = file.split('.pkl')[0]
            if color_by_conservation == 'all' or color_by_conservation == 'centers':
                identify_conserved_water_clusters(networks, centers, dist_cutoff=1.0, filename_base=f'{cluster_filebase}_{name}_conservation')
            if color_by_conservation == 'all' or color_by_conservation == 'connections':
                identify_conserved_water_interactions_clustering(networks, centers, max_connection_distance=3.0, dist_cutoff=1.0, filename_base=f'{cluster_filebase}_{name}_conservation')

    if classify_waters:
        csvs = [f for f in os.listdir(csv_dir) if f.endswith('.csv')]
        residue_analysis.plot_interactions_from_angles(csvs)


 
def parse_arguments():
    """
    Parse arguments from command line
    """
    parser = argparse.ArgumentParser(description='Perform analysis using WatCon')
    parser.add_argument('--input', type=str, help='Input file', default=None)
    parser.add_argument('--analysis', type=str, help="Analysis Input File", default=None)
    parser.add_argument('--name', type=str, help='Identifiable name', default='results')
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()
    if args.input is not None:
        structure_type, kwargs = parse_inputs(args.input)
        results = run_watcon(structure_type, kwargs)
        os.makedirs('watcon_output', exist_ok=True)
        with open(f'watcon_output/{args.name}.pkl', 'wb') as f:
            pickle.dump(results, f)
    
    if args.analysis is not None:
        kwargs = parse_analysis(args.analysis)
        run_watcon_postanalysis(**kwargs)


    
