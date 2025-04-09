Frequently Asked Questions
==========================


.. dropdown:: Do I need a modeller license to use WatCon?

   Core WatCon functionality relies on the Modeller package. Although Modeller is not required for all WatCon analysis, the :mod:`WatCon.sequence_processing` module requires a working Modeller installation. Details on installing Modeller (and obtaining the license key) can be found `here <https://salilab.org/modeller/>`_ . 

.. dropdown:: How do I align my structures if I have multimers?

   .. role:: python(code)
       :language: python

   WatCon interfaces with Modeller in order to align structures based on sequence. This is done using the :func:`WatCon.sequence_processing.perform_structure_alignment` function. By default, alignmemts will only use the sequence from Chain 'A', but can be modified by the user with the :python:`same_chain` flag. Then, once the chains of interest have been aligned, the :func:`WatCon.sequence_processing.align_with_waters` function can be used to transform all corresponding water/ions/ligands/chains to the new coordinate space. By default, only the chains of interest are going to be aligned, but the user can change the :python:`selected_chain_only` flag to ensure that all atoms are translated and rotated according to the new coordinate space. 

.. dropdown:: Why am I getting conservation scores of 0 for all of my structures?

   This is likely due to a misalignment with your clusters and your structures. Refer to the :doc:`Combining Data <combining_different_data.rst>` section for details on how to do this correctly. 

.. dropdown:: How can I transfer my cluster coordinates onto a different set of related structures?

   The :mod:`WatCon.sequence_processing` module is not only useful for PDB processing. The :func:`WatCon.sequence_processing.perform_structure_alignment` function is also useful for moving between different coordinate alignment spaces. A list can be inputted into this function (rather than a directory), and all structures from that list will be aligned. The transformation matrices returned from this function can be used to move from one coordinate space to another. For instance, assume that a set of clusters were created from one list of structures (set A) and you want to compute a conservation score for each PDB in set B to this set of clusters in set A. This can be accomplished by transforming the coordinates of the clusters in set A to match the coordinates in set B.

   .. code-block:: python
      
      from WatCon import sequence_processing, find_conserved_networks

      #Cluster file of interest
      cluster_file = 'CLUSTER_setA.pdb'
      cluster_coordiantes = find_conserved_networks.get_coordinates_from_topology(cluster_file)

      #Set of PDBs to align. All structures will be aligned to the FIRST structure in the list
      pdbs = ['setA.pdb', 'setB.pdb']

      #Calculate rotation information. IMPORTANT: Set sort_pdbs=False in order to ensure that setA.pdb is aligned to setB.pdb
      rotation_information = sequence_processing.perform_structural_alignment(pdbs, out_dir='aligned_pdbs_new', sort_pdbs=False)

      #Transform cluster coordinates
      transformed_centers = np.dot(cluster_coordinates, rotation_information['Rot'][0].T) + rotation_information['Trans'][0]


.. dropdown:: Why are my fasta files not being recognized?

   WatCon relies on a standard naming scheme for recognizing the correct sequence in an alignment file and the correct corresponding fasta file. Refer to the :doc:`Organizing files section <organizing_files>` for more details on the WatCon-anticipated naming schemes of files.

.. dropdown:: Why are my atom number selections off by 1 index?

   WatCon uses a 0-based atom indexing scheme, while most structures use a 1-based indexing scheme. However, outputted WatCon PyMOL files are created with this offset in mind. Therefore, you only need to adjust atom numbers by 1 if take atom indexes directly from WatCon output files.

.. dropdown:: How do I know what clustering method to use?

   Choosing the most appropriate clustering method is an involved problem that is not always trivial to answer. We have achieved the most accurate clustering results from using the HDBSCAN methods, as this method requires the fewest number of parameters and therefore is more robust to changes in water distributions. However, HDBSCAN can become quite slow, and so OPTICS or DBSCAN can be preferable for systems where you know approximately how spaced apart your centroids should be. We recommend, especially for crystal structures, plotting the cluster centers over all of the aligned structure files to see if the locations of the clustered waters are consistent with water clusters that can be visually observed from the aligned structures.

.. dropdown:: Why can't I see my densities?

   The :func:`WatCon.generate_dynamic_networks.collect_densities` function will output a .dx file which can be visualized with PyMOL and other visualization softwares. In PyMOL, load the .dx file and click the :python:`A->mesh/surface/volume` options in order to see the densities. If densities are not visible, then it is likely that your histograms were too sparse. If this is the case, check your simulation files to make sure that you have a sufficient number of frames (you may need multiple hundreds or thousands of frames to achieve appropriate sampling) and that your selected water region is correct.


.. dropdown:: Do I have to use input files to run WatCon analysis?
   
   No! If you prefer, you can use the :doc:`Python API <../api>` directly.

.. dropdown:: Why am I getting OOM errors?

   There is a chance that, especially when working with large trajectories, WatCon analysis will use too much working memory and crash. This should only happen if you set the :python:`return_network=True` flag in the :func:`WatCon.generate_dynamic_networks.initialize_network` function with a sufficiently large active region selected. By default, this flag is set to False to avoid these issues.

.. dropdown:: How do I know what two references to use for my angle calculations?

   The two angle references are generally arbitrary, and can be any two points in 3D space. Generally, you want these points to be noticeably far away from each other (broadly speaking, if the references are too close then the two angles will be calculated as too similar and a useful distribution of angles will not be achieved. If processing multiple sets of unaligned structures, then we recommend using two reference atoms which are generally very static, as these references can then be consistent across multiple sets of structures, even if the coordinate space is different.

.. dropdown:: Why are the interactions on my PyMOL projections clearly inaccurate?

   This can happen under a number of circumstances, but normally occurs when structures loaded are slightly different than those inputted into WatCon initially. For trajectories, ensure that the .pml file corresponds to the frame that you are currently viewing on PyMOL. Additionally, a sufficiently large active region could potentially cause interactions to be found which cross the periodic boundaries of the simulation box. If this happens, the resulting connections will look incorrect, but are consistent with periodic boundary condition definitions. Finally, ensure that the atom labeling for the structure file loaded in PyMOL is the same as that inputted into WatCon, differences in input file structure (even if the underlying protein is the same) can cause changes in atom numbers, which will affect projection accuracy. 


