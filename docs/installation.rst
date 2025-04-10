Installation Guide
==================

Instructions to install WatCon are as follows:


1. Clone the WatCon Repository
------------------------------

Clone the WatCon repository and change your working directory into WatCon

.. code-block:: bash

   git clone https://github.com/kamerlinlab/WatCon.git

   cd WatCOn


2. Create a WatCon Conda Environment
------------------------------------

.. code-block:: bash

   conda env create -f WatCon.yaml


Once completed, a WatCon environment should have been created so that you can freely use the code. With input files, call WatCon on the command line by

.. code-block:: console

   $ python -m WatCon.WatCon --input input.txt --name name_of_system


Note
----

Certain features of WatCon are dependent on the Modeller package by the Sali Lab. This package requires a license which can be obtained `here <https://salilab.org/modeller/>`_. 
