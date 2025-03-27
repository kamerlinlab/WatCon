Combining Different Sources of Data
-----------------------------------


There are multiple ways at combining data from different sources for WatCon analysis. We provide a series of examples to demonstrate the recommended methods.


Combining Multiple Replicas from Simulation Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you are combining results from multiple replicas of the same protein, we recommend to concatenate all desired trajectories (taking care to align to a common reference frame) as part of your preprocessing procedure before running WatCon. 


Combining Different Collections of Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Since WatCon is divided into separate dynamic and static processing modules, care needs to be taken to combine static and dynamic data. Similarly, if two sets of static data are run separately and combined post-WatCon analysis, care needs to be taken to ensure proper alignment. 

TALK ABOUT CLUSTERING AND ROTATION MATRICES 
