Organizing Files
----------------


WatCon operates under the assumption that the user has a common directory structure. 
Although naming schemes are generally flexible to the users' preference, we provide an example structure that can be utilized.


Directory Structure
~~~~~~~~~~~~~~~~~~~

The recommended initial structure for a static structure analysis is as follows:

.. code-block:: txt

    WatCon_Analysis_Folder
    ├── pdbs
    ├── clean_pdbs #Files created using pdb4amber, see :doc:`Basic Tutorial <basic_tutorial>`
    ├── fasta
    alignment.txt #MSA alignment file
    all_fastas.fa #Combined fasta file
    input_file.txt #Input file (optional)

After using WatCon to align structures, the resulting directory structure should look like this:

.. code-block:: txt

    WatCon_Analysis_Folder
    ├── pdbs
    ├── clean_pdbs
    ├── aligned_with_waters
    ├── fasta
    alignment.txt
    all_fastas.fa
    input_file.txt

After WatCon analysis, the resulting directory structure will look like this:

.. code-block:: txt

    WatCon_Analysis_Folder
    ├── pdbs
    ├── clean_pdbs
    ├── aligned_with_waters
    ├── fasta
    ├── watcon_output #If using input files
    ├── cluster_pdbs #If clustering
    ├── msa_classification #If using water-angle classification
    ├── pymol_projections #If creating PyMOL projections
    alignment.txt
    all_fastas.fa
    input_file.txt


Naming Schemes
~~~~~~~~~~~~~~

WatCon assumes a consistent naming scheme between structures, trajectories, and fasta files. 

WatCon assumes that aligned PDB files will follow the naming convention of ${PDB_ID}_aligned.pdb and that fasta files follow the naming convention of ${PDB_ID}.fa. 
