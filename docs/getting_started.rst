Getting Started
===============


Running WatCon
--------------

WatCon has been designed to be implemented via the use of either an input file or direct python interface.


Sample Input File
~~~~~~~~~~~~~~~~

An example of a WatCon input file is provided below. 

.. code-block:: txt
    
   ; WatCon Sample Input File
   
   ; Initialization
   structure_type: dynamic                          ; Create static or dynamic networks
                                                    ;     Use dynamic for multiple pdb models
   structure_directory: structure_dir               ; Directory of structure files 
   topology_file: WT_PTP1B_Unlig_Open.gro             ; Name of topology (required for dynamic)
   trajectory_file: run_1.xtc                       ; Name of trajectory (required for dynamic)
   network_type: water-protein                      ; Use only waters (water-water) or also protein atoms (water-protein)
   include_hydrogens: off                           ; If include_hydrogens = on, create a directed graph
                                                    ;     MDAnalysis 'protein' selection
   water_name: default                              ; Any custom water names
   multi_model_pdb: False                           ; pdb files have multiple models (typical of NMR structures)
   max_distance: 3.0                                ; Max distance between two atoms to be considered
                                                    ;     in an HBond (recommended 3.3 if using static 
                                                    ;     structures with no hydrogens, 1.8 if dynamic   
                                                    ;     structures with hydrogens)
   angle_criteria: 150                              ; Specify criteria for calculating HBonds with 
                                                    ;     angles+distances (recommended 120 if hydrogens
                                                    ;     are present)
   
   ; Property calculation
   density: on
   connected_components: on
   interaction_counts: on
   per_residue_interactions: on
   characteristic_path_length: on
   graph_entropy: on
   clustering_coefficient: on
   save_coordinates: on
   analysis_selection: all                          ; Selection for analysis 
                                                    ;     (all, active_site, not_active_site)
   
   ; Active site definition
   active_site_reference: resid#220#and#name#CA     ; MDAnalysis selection language to center active site
   
   active_site_only: on                             ; Indicate whether to only calculate water networks
                                                    ;     around an active site 
   active_site_radius: 11                           ; Radius of active site around refernce
   
   ; Visualization
   project_networks: off                            ; Create PyMOL files per pdb/frame
   
   ; Clustering
   cluster_coordinates: on                          ; Perform a clustering analysis on all coordiantes
   clustering_method: hdbscan                       ; Clustering algorithm: dbscan, hdbscan, or optics
   min_cluster_samples: 100                         ; Minimum samples per cluster
   eps: 0.0                                         ; Eps value for clustering
   
   ; MSA Indexing
   msa_indexing: on                                 ; Utilize/perform an MSA
   alignment_file: alignment.txt                    ; Name of alignment file (if file does not exist, 
                                                    ;     Modeller will be used to write this file)
   make_fastas: off                                 ; If 'on', WatCon will make fasta files 
                                                    ;     from the pdbs in structure_directory
   combined_fasta: all_fastas.fa                    ; Name of combined fasta file
   fasta_directory: fasta                           ; Directory containing individual fasta files
   MSA_reference_pdb: WT_PTP1B_Unlig_Open.gro         ; Any pdb which can be used as a reference
                                                    ;     (active_site_reference needs to be accurate 
                                                    ;     for this structure)
   
   ; Classify waters from MSA
   classify_water: on                               ; Classify water by angles and MSA
   water_reference_resids: 70,#153                  ; Residue positions to use as reference points 
                                                    ;     (in relation to MSA_reference_pdb)
   classification_file_base: DYNAMIC_OPEN_1         ; Classification file basename
   
   ; Miscellaneous
   num_workers: 8                                   ; Number of cores available for parallelization


Execution
~~~~~~~~~

WatCon can then be executed by the command

.. code-block:: console

   $ python -m WatCon.WatCon --input input_file.txt --name name_of_system


Which will output any PDB files and PyMOL files as specified by the user. Results will be outputted in a .pkl file which can then be loaded and analyzed further.


Analyzing calculated metrics
----------------------------



Visualizing Results with PyMOL
------------------------------




