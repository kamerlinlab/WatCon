Combining Different Sources of Data
-----------------------------------


There are multiple ways at combining data from different sources for WatCon analysis. Examples on instances of combination are shown in the :docs:`Tutorials <../tutorials>`.


Combining Multiple Replicas from Simulation Data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you are combining results from multiple replicas of the same protein, we recommend using WatCon in one of two ways:

1. Concatenate all desired trajectories (taking care to align to a common reference frame) as part of your preprocessing procedure before running WatCon. 
2. Run WatCon individually for each trajectory and combine results after.


Combining Different Collections of Data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Since WatCon is divided into separate dynamic and static processing modules, care needs to be taken to appropriately combine static and dynamic data. Similarly, if two sets of static data are run separately and combined post-WatCon analysis, care needs to be taken to ensure proper alignment. Generally the :mod:`WatCon.sequence_processing.perform_structure_alignment` function is very useful in manipulating coordinates among related but not aligned structures. We recommend the general procedure for switching between coordinate spaces:

.. code-block:: python

   from WatCon import sequence_processing

   coords = #Coordinates from clusters, waters, protein_atoms, etc. (aligned to the second pdb in the list below)
   pdbs = ['target.pdb', 'source.pdb']

   #Get rotation information from alignment
   rotation_information = sequence_processing.perform_structure_alignment(pdbs, sort_pdbs=False)

   #Transform coordinates
   transformed_coords = np.dot(coords, rotation_information['Rot'][0].T) + rotation_information['Trans'][0])


WatCon will output lists of dictionaries of calculated values. If using :mod:`WatCon.generate_dynamic_networks`, you will receive a tuple of lists all with one item per trajectory frame. If using :mod:`WatCon.generate_static_networks` instead, you will receive a tuple of lists with one item per structure file. 
