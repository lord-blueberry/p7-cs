# p7-cs
Repository contains the work done in the P7 Project, "Compressed Sensing Image Reconstruction for CASA". The Folder ./casa-cs/casa_impl  contains the source code which can be run in CASA.

The code was integrated in CASA 5.3.0. Older versions may run, but it was not tested.

## Installation Guide
CASA_ROOT 
* Copy file ./casa-cs/casa_impl/p7_cs.py into  CASA_ROOT/lib/python2.7
* apply patch on CASA_ROOT/lib/python2.7/tasks.py. The patch file is in ./casa-cs/casa_impl/tasks.py.patch
* Download Gurobi
* make sure you have a Gurobi licencse
* install Gurobi python bindings to CASA. Call CASA_ROOT/bin/python -m GUROBI_ROOT/linux64/setup.py
* Add GUROBI_ROOT/linux64/bin to CASALD_LIBRARY_PATH environment variable. (modify your .bashrc file for ease of access)

* Upgrade Numpy version of CASA:
  - call CASA_ROOT/bin/python -m pip install --upgrade pip setuptools wheel
  - call CASA_ROOT/bin/python -m pip install numpy --upgrad









