## What is this github project about?
This github site describes various modified versions of the GFDL Idealized Spectral Atmospheric Models developed in [Prof Gang Chen's research group at UCLA](http://gchenpu.com/group/).  These modifications are designed to answer certain research questions in atmospheric and climate dynamics that cannot be easily addressed by playing with the default model code.  As different research groups may have different versions of the idealized models, we have provided the comparison between default and modified configurations to facilitate the use of these modifications in other model versions or different models.  The models are also run through a python interface written and tested at the NCAR supercomputer.  Questions about these modifications can be sent to *gchenpu at atmos dot ucla dot edu*.

## Where to start:
Documentation of the GFDL Idealized Global Atmospheric Models with Spectral Dynamics can be found here:
http://www.gfdl.noaa.gov/idealized-spectral-models-quickstart

## Specific configurations
- [how to test the spectral transforms in the model?](../exp/f90/spectral_tutorial)
- [how to set up a zonally symmetric version of the model with the same climatology of the full model?](../exp/hs/zonal_symmetric_model)
- more to be added

## Summary of directory contents
bin/  
- mkmf script for creating Makefiles  
- template files for the mkmf, compile, and run scripts  
- time_stamp.csh used by run script  

docs/: documentations  

exps/  
- different configurations of the models
  - f90: Examples of very simple models using the subroutine in the spectral model
  - hs: Held Suarez model

fmspy/: a python interface that helps compile and run different models

netcdf/: netcdf libraries used

postprocessing/  
- mppnccombine: a tool for combining distributed diagnostic output files
- plevel: a tool for interpolating data from model levels to pressure levels

src/: default source code files

run_fms.py:  python script to compile and run the model
