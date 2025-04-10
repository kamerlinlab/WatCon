Getting Started
===============


Preparing Structures
--------------------

WatCon is an analysis code created for the purposes of water network detection, analysis, and comparison across static and dynamic structures. Directions on how to prepare structures for analysis are as follows.


Static Structures
~~~~~~~~~~~~~~~~~

We require that all static structures be aligned before WatCon analysis. We provide the functionality to perform this alignment using the :mod:`WatCon.sequence_processing` module, which requires structures in PDB format with no headers. This formatting style can be easily accomplished using the ``pdb4amber`` function in the `AmberTools suite <https://ambermd.org/AmberTools.php>`_. 


Dynamic Structures
~~~~~~~~~~~~~~~~~~

We require all dynamic (trajectories) to have been post-processed, with any PBC errors resolved and protein coordinates aligned. Since there are numerous methods at accomplishing this, we leave this process up to the user. To read dynamic information, WatCon leverages the use of `MDAnalysis <https://www.mdanalysis.org/>`_, and so can read any MDAnalysis-readable topologies and trajectories. 


Multiple Sequence Alignment
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Numerous WatCon functionality is dependent on a common sequence indexing across all structures. We require a multiple sequence alignment (MSA) to be conducted and written to a file in PIR format for this purpose. We have a basic sequence alignment tool implemented in the :mod:`WatCon.sequence_processing` module, which leverages `Modeller <https://salilab.org/modeller/>`_, but we recommend performing a more rigorous alignment for structures which are not closely related and always recommend visually inspecting alignment outputs to ensure accuracy.



Running WatCon
--------------

WatCon has been designed to be implemented via the use of either an input file or direct python interface.


Sample Input File
~~~~~~~~~~~~~~~~~

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
   include_hydrogens: on                            ; If include_hydrogens = on, create a directed graph
                                                    ;     MDAnalysis 'protein' selection
   water_name: default                              ; Any custom water names
   multi_model_pdb: False                           ; pdb files have multiple models (typical of NMR structures)
   max_distance: 3.0                                ; Max distance between two atoms to be considered
                                                    ;     in an HBond (recommended 3.3 if using static 
                                                    ;     structures with no hydrogens, 1.8 if dynamic   
                                                    ;     structures with hydrogens)
   angle_criteria: 150                              ; Specify criteria for calculating HBonds with 
                                                    ;     angles+distances (recommended 120 if hydrogens
                                                    ;     are present -- only activates if include_hydrogens=on)
   
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
                                                    ;     (all, active_region, not_active_region)
   
   ; Active site definition
   active_region_reference: resid#220#and#name#CA     ; MDAnalysis selection language to center active site
   
   active_region_only: on                             ; Indicate whether to only calculate water networks
                                                    ;     around an active site 
   active_region_radius: 11                           ; Radius of active site around refernce
   
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
   combined_fasta: all_fastas.fa                    ; Name of combined fasta file
   fasta_directory: fasta                           ; Directory containing individual fasta files
   MSA_reference_pdb: WT_PTP1B_Unlig_Open.gro         ; Any pdb which can be used as a reference
                                                    ;     (active_region_reference needs to be accurate 
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


Analyzing Calculated Metrics
----------------------------

To increase ease in combining results across multiple trajectories or multiple static structures, we also allow for supplemental analysis to be conducted via direct python interface or input files following this construction:


.. code-block:: txt

    ; Sample input file for WatCon analysis
    
    
    ; Initialize
    concatenate: PTP1B_closed_1,PTP1B_closed_2,... ; Indicate which runs to concatenate. All others will be treated separately
    input_directory: watcon_output                 ; Folder which contains outputted WatCon .pkl files
    
    
    ; Basic metric analysis
    histogram_metrics: on                          ; Will make basic matplotlib histograms of metrics according to desired concatenation                
    residue_interactions: on                       ; Will make bar graphs of per-residue water interaction frequency
    reference_topology: WT_PTP1B_Closed.gro        ; Reference topology for residue_interaction analysis
    interaction_cutoff: 0.0                        ; Cutoff of interaction score for writing residue_interaction bar graph    
    
    ; Density analysis
    calculate_densities: on                        ; Use MDAnalysis density to output density of waters (ONLY FOR DYNAMIC_NETWORKS)    
    active_region_definition: around 11 (resid 220 and name CA) ; MDAnalysis selection language for active site
    
    ; Cluster conservation
    cluster_filebase: STATIC_CLOSED                     ; Name of pdb file containing clusters to compare to
    cluster_concatenated: off                      ; Indicate whether to cluster concatenated coordinates
    calculate_commonality: bar                     ; Produce commonality plot either as a 'bar' bar graph or 'hist' histogram 
    color_by_conservation: all                     ; Produce .pml file coloring either 'centers', 'connections' or 'all' by conservation
    
    ; Residue-water classification
    classify_waters: on                            ; Use outputted .csv files from 2-angle classification to generate scatter/density plots

WatCon post-analysis can then be executed by the command

.. code-block:: console

   $ python -m WatCon.WatCon --analysis analysis_file.txt


We further note that both an input file and analysis file can be passed simulataneously, i.e:

.. code-block:: console

   $ python -m WatCon.WatCon --input input_file.txt --analysis analysis_file.txt --name name_of_system


Visualizing Results with PyMOL
------------------------------

Numerous modules will output files to visualize water networks or clustered water positions. We provide a brief description of common file types and how to visualize in PyMOL

* Cluster centers (.pdb): Cluster centers will often be saved as PDB files with dummy water oxygen atoms. These can safely be loaded in PyMOL similar to any other PDB file
* Snapshot water networks (.pml): PyMOL (.pml) files containing information regarding connections among waters are often written, and can easily be loaded in PyMOL by first loading the corresponding structure file and then typing :code:`@FILE.pml` in the PyMOL command line

  .. image:: images/open_projection.png
     :width: 400
     :align: center

  .. note:: When loading .pml files for trajectory frames, we recommend only loading in the trajectory frame of interest into PyMOL first, and then follow by loading the .pml. This will reduce unecessary wait time when loading the connections.

* Density distributions (.dx): Density distributions can be loaded into PyMOL with the :code:`load` command, similar to PDB files. For most clear results, load a corresponding structure file first before loading the densities. 

  .. note:: When loading .dx files calculated from MDAnalysis, a structure file written from that same script should be loaded first, to guarantee proper alignment


