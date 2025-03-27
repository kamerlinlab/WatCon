Using the Python API versus Input Files
---------------------------------------

We allow the user to run WatCon using input files or with the Python API directly. When using input files, an outputted ".pkl" file is generated which can be reloaded to perform subsequent analysis. Although WatCon analysis is very fast for static structures, it can take several minutes to process a large trajectory, and so the outputted .pkl files are very useful for data analysis as WatCon does not need to be rerun every time a comparison is being made. Simulataneously, the Python API can be more convenient to use, as the outputted data can be accessed and manipulated directly. 


When utilizing input files (example given in the :doc:`Getting Started <../getting_started>`) section, results can be loaded and manipulated as follows:

.. code-block:: python

   import pickle

   with open('WATCON_OUTPUT.pkl', 'rb') as FILE:
       data = pickle.load(FILE)


   metrics_dict, networks, cluster_centers, names = data


And the Python API can be utilized in the same manner:

.. code-block:: python
   
   import WatCon.generate_static_networks as generate_static_networks
 
   metrics_dict, networks, cluster_centers, names = generate_static_networks.initialize_network(\*\*kwargs)


Similarly, post analysis can be run with input files or using the Python API directly. Post analysis calculated using input files are all written to disk as image files or pdb/pml files for molecular visualization. Utilizing the Python API can provide greater flexibility in post analysis.


 
