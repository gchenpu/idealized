## Zonally symmetric model modified from the [default Held-Suarez model](../default)

## modified changes to the code
- tendency option to save the tendencies from integrating the model by one time step; 
- forcing option to add the prescribed forcing to the model

## see the changes in the code
- default src code: [spectral_dynamics.F90](../../../src/atmos_spectral/model/spectral_dynamics.F90)
- modified src code: [spectral_dynamics.F90](./src/spectral_dynamics.F90)

## modified namelist parameters
- set the model truncation to zonally symmetric (i.e., fourier_inc = num_fourier)
- create climatology.nc using matlab; read the file as the initial condition (i.e., initial_state_option    = 'input'); set up the option in ic_from_external_file_nml for the input file

## Summary of directory contents
run_ctl_forcing/  
- model configuration to obtain the eddy forcing

run_ctl/  
- model configuration with prescribed eddy forcing

src/  
- pathnames of the source code 
