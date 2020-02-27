# -----  zonally symmetric model modified from the Held-Suarez model ----

updates in the code
diff ../../../src/atmos_spectral/model/spectral_dynamics.F90 ./src/spectral_dynamics.F90

steps:
1. set the model truncation to zonally symmetric (i.e., fourier_inc = num_fourier)

2. create climatology.nc using matlab; 
   read the file as the initial condition (i.e., initial_state_option    = 'input'); 
   setup the option in ic_from_external_file_nml for the input file

3. tendency option to save the tendencies from integrating the model by one time step; 
   forcing option to add the prescribed forcing to the model
