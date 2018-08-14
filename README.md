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

* Upgrade Numpy version of CASA (not needed, but old numpy versions can increase the runtime a lot):
  - call CASA_ROOT/bin/python -m pip install --upgrade pip setuptools wheel
  - call CASA_ROOT/bin/python -m pip install numpy --upgrade
  

## Usage
In the CASA IPython interface (call CASA_ROOT/bin/casa), you can use the new p7_cs task. The usage is similar to the tclean task, albeit with additional parameters. Also, since it is a proof of concept implementation, not all ways in which tclean can be called are supported.

For an example dataset, download [SNR_G55_10s.calib.ms from here](https://casaguides.nrao.edu/index.php/VLA_CASA_Imaging-CASA5.0.0).

### Additional parameters compared to tclean
*cs_alg= String*, this specifies the algorithm. Valid input:
* positive_deconv
* L1
* L2
* L1+L2
* haar
* TV
* starlet
(defaults to starlet if none specified)

*lambda_cs=Float*, Specifies the lambda value for Compressed Sensing. If none is given, defaults to 0.05

*psf_threshold=Float* used to truncate small parts of the PSF. Valid range is betweeen 0 and 1

*psf_cutoff=False* used to reduce the psf to half it's size. This was not used in the final run

*starlet_levels=Integer* How many stalets are used. Starlets has a multi-resolution representation. The smallest starlet (Level 1) is 5 pixels wide. Each more level doubles the starlet size (Level 2 is 10 pixels, level 3 starlet is 20 pixels wide). Significantly impacts the runtime, level 3 was the maximum used in this project.

*lambda_estimate=["XXX__objectiveVal.csv", "XXX_deconv_solution.csv"]* used for the miller lambda estimation. The algorithm "positive_deconv" produces two extra csv files: "XXX__objectiveVal.csv" and "XXX_deconv_solution.csv". They are needed for the miller lambda estimation.

### Example Usage
deconvolve image and create basis for miller lambda estimation

p7_cs(vis='SNR_G55_10s.calib.ms',niter=0,imagename="deconv", imsize=[128,128], cell='26arcsec', interactive=False, stokes='I', cs_alg='positive_deconv', psf_threshold=0.02, psf_cutoff=False, wprojplanes = 1, pblimit=-1.0)

starlet regularization used in the project

p7_cs(vis='SNR_G55_10s.calib.ms',niter=0,imagename='G55_starlets', imsize=[128,128], cell='26arcsec', interactive=False, stokes='I', cs_alg='starlets', starlet_levels=3, psf_threshold=0.02, psf_cutoff=False, lambda_estimate=["deconv_objectiveVal.csv", "deconv_solution.csv"], wprojplanes = 1, pblimit=-1.0)









