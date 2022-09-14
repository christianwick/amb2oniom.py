# Amber2Oniom.py // amb2oniom.py 
------------------------------

Convert Amber topologies to oniom input files.


## Description:
------------

amb2oniom.py is a python3 script which converts amber topologies to oniom input
files, which include all parameters from the amber topolog rather then
relying on the generic parameter sets avaiable in gaussian. 
All the amber topology functionality is added via ParmEd
(https://github.com/ParmEd/ParmEd) and a working parmed installation is
necessary to run amb2oniom.py. It is recommended to use the latest parmed
version (or at least a stable version past 3.0.0).


## Install amb2oniom.py:
---------------------

amb2oniom.py can be copied to any directory and used out of the box as 
long as you have a running python3 installation and parmed avaiable (either
from AmberTools: http://ambermd.org or directly from github:
https://github.com/ParmEd/ParmEd )


## What it can do:
---------------

  - write gaussian input files for oniom calculations with Amber force fields.
  - assign hi, medium and low layers based on ambmask like strings
  - strip the topology and write coordinates and toplologies
  - detect link atoms (atom types need to be assigned by hand)
  - select active atoms (the belly region)
  
  
## What you need:
--------------

 - a correct amber topology and coordinate files
 - a working copy of parmed (it is e.g. included in the AmberTools)
 

## Examples:
---------

  - print help message and exit
    
    ./amb2oniom.py -h

  - read in amb.top and amb.crd, select residue :20 as hi layer and 
    set the belly region to 10A around the high layer:

    ./amb2oniom.py -p amb.top -c amb.top -hi ":20" -act ":20<:10"

  - read in amb.top and amb.ncrst, strip resdiues 100-150:

    ./amb2oniom.py -p amb.top -c amb.rst7 -strip ":100-150"

  - modify gaussian route section and change charge multiplity

    ./amb2oniom.py -p amb.top -c amb.crd -route "#P ONIOM(AM1/AMBER=Softonly)"\
          -cm "0 1 0 2 0 2"

NOTE: Don't forget to add the link atom types to your link atoms!! 
      amb2oniom.py prints atom types of all atoms connected to the 
      link atom host as a guide.
 
 
## License / Permission of use
The code is provided as-is without any warranty for its correctness under the MIT license (see LICENSE file)

## Contact
Dr. Christian R. Wick
Friedrich-Alexander-Universität Erlangen-Nürnberg,  PULS Group, Institute for theoretical Physics and Competence Unit for Scientific Computing (CSC), Interdisciplinary Center for Nanostructured Films (IZNF), Cauerstrasse 3, 91058 Erlangen, Germany
 
 
 
 
