# p7-cs
Repository contains the work done in the P7 Project, "Compressed Sensing Image Reconstruction for CASA". The Folder ./casa-cs/casa_impl  contains the source code which can be run in CASA.

The code was integrated in CASA 5.3.0. Older versions may run, but it was not tested.

## Installation Guide
CASA_ROOT is the installation directory of casa. If you [downloaded CASA 5.3.0](https://casa.nrao.edu/casa_obtaining.shtml), the folder is called casa-release-5.3.0-XXX.elX
* Copy file ./casa-cs/casa_impl/p7_cs.py into  CASA_ROOT/lib/python2.7
* apply patch on CASA_ROOT/lib/python2.7/tasks.py. The patch file is in ./casa-cs/casa_impl/tasks.py.patch
* Download Gurobi
* make sure you have a Gurobi licencse
* install Gurobi python bindings to CASA. Call CASA_ROOT/bin/python -m GUROBI_ROOT/linux64/setup.py
* Add GUROBI_ROOT/linux64/bin to CASALD_LIBRARY_PATH environment variable. (modify your .bashrc file for ease of access)

* Upgrade Numpy version of CASA:
  - call CASA_ROOT/bin/python -m pip install --upgrade pip setuptools wheel
  - call CASA_ROOT/bin/python -m pip install numpy --upgrade
  

## Usage
In the CASA IPython interface (call CASA_ROOT/bin/casa), you can use the new p7_cs task. The usage is similar to the tclean task, albeit with additional parameters. Also, since it is a proof of concept implementation, not all ways in which tclean can be called are supported.

For an example dataset, download [SNR_G55_10s.calib.ms from here](https://casaguides.nrao.edu/index.php/VLA_CASA_Imaging-CASA5.0.0).

### Additional parameters compared to tclean
cs_alg=None

lambda_cs=None

psf_threshold=None

psf_cutoff=False

starlet_levels=Integer



lambda_estimate=None

### Example Usage
p7_cs(vis='SNR_G55_10s.calib.ms',niter=0,imagename='G55_starlets', imsize=[128,128], cell='26arcsec', interactive=False, stokes='I', cs_alg='starlets', psf_threshold=0.02, psf_cutoff=False, lambda_estimate=lamb_est, wprojplanes = 1, pblimit=-1.0, starlet_levels=3)









