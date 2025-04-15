Installation Guide
==================

Instructions to install WatCon are as follows:


1. Clone the WatCon Repository
------------------------------

Clone the WatCon repository and change your working directory into WatCon

.. code-block:: bash

   git clone https://github.com/kamerlinlab/WatCon.git

   cd WatCon


2. Create a WatCon Conda Environment
------------------------------------

.. code-block:: bash

   conda env create -f WatCon.yaml


3. Add Modeller license
-----------------------

After you create the WatCon conda environment, you will receive this message:

.. code-block:: txt

   Edit /anaconda3/envs/WatCon/lib/modeller-10.6/modlib/modeller/config.py
   and replace XXXX with your Modeller license key
   (or set the KEY_MODELLER environment variable before running 'conda install').

This message is prompted because certain features of WatCon are dependent on the Modeller package by the Sali Lab. This package requires a license which can be obtained `here <https://salilab.org/modeller/>`_. Once you receive the license, simply edit the modeller config file with the license key.


Once completed, a WatCon environment should have been created so that you can freely use the code. With input files, call WatCon on the command line by

.. code-block:: console

   $ python -m WatCon.WatCon --input input.txt --name name_of_system

or simply import WatCon as a python package directly:

.. code-block:: python

   import WatCon.sequence_processing

   sequence_processing.pdb_to_fastas('clean_pdbs/structure1.pdb', 'fasta', 'PTP1B')
