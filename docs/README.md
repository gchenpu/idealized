Where to start
~~~~~~~~~~~~~~
Documentation of these models can be found here:
http://www.gfdl.noaa.gov/idealized-spectral-models-quickstart

Summary of directory contents
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bin/
     -- mkmf script for creating Makefiles
     -- template files for the mkmf script
     -- list_paths generates a list of source code files using find
     -- time_stamp.csh used by run script

barotropic/ -- compile script. (uses mkmf)
            -- list of source code files needed for barotropic model. (used by mkmf)
            -- script for running the barotropic model.
            -- input data files needed for barotropic model.

shallow/   -- compile script. (uses mkmf)
           -- list of source code files needed for shallow water model. (used by mkmf)
           -- script for running the shallow water model.
           -- input data files needed for shallow water model.

full_dynamics/ -- compile script. (uses mkmf)
               -- list of source code files needed for full spectral dynamics model. (used by mkmf)
               -- script for running the full spectral dynamics model.
               -- input data files needed for full spectral dynamics model.

src/ -- source code files

postprocessing/ -- contains source code for a tool for combining distributed diagnostic output files
