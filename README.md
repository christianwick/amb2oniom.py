Amber2Oniom.py // amb2oniom.py 
------------------------------

Convert Amber topologies to oniom input files.


Description:
------------

amb2oniom.py is a python3 script which converts amber topologies to oniom input
files, which include all parameters from the amber topolog rather then
relying on the generic parameter sets avaiable in gaussian. 
All the amber topology functionality is added via ParmEd
(https://github.com/ParmEd/ParmEd) and a working parmed installation is
necessary to run amb2oniom.py. It is recommended to use the latest parmed
version (or at least a stable version past 3.0.0).


Install amb2oniom.py:
---------------------

amb2oniom.py can be copied to any directory and used out of the box as 
long as you have a running python3 installation and parmed avaiable (either
from AmberTools: http://ambermd.org or directly from github:
https://github.com/ParmEd/ParmEd )


What it can do:
---------------

  - write gaussian input files for oniom calculations with Amber force fields.
  - assign hi, medium and low layers based on ambmask like strings
  - strip the topology and write coordinates and toplologies
  - detect link atoms (atom types need to be assigned by hand)
  - select active atoms (the belly region)
  
  
What you need:
--------------

 - a correct amber topology and coordinate files
 - a working copy of parmed (it is e.g. included in the AmberTools)
 
 
 
Contributors:
-------------
 
  - Chris Wick (principle developer)
 
 
 
 