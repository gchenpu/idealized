## Zonally symmetric model modified from the [default Held-Suarez model](../default)  
The method follows the appendix of  
Kushner, P. J., & Polvani, L. M. L. (2004). Stratosphere–Troposphere Coupling in a Relatively Simple AGCM: The Role of Eddies. Journal of Climate, 17(3), 629–639. [Link](https://doi.org/10.1175/1520-0442(2004)017<0629:SCIARS>2.0.CO;2)  

An extension of this method was decribed in  
Domeisen, D. I. V., Sun, L., & Chen, G. (2013). The role of synoptic eddies in the tropospheric response to stratospheric variability. Geophysical Research Letters, 40(18), 4933–4937. [link](https://doi.org/10.1002/grl.50943)

Modified changes to the code
- tendency option to save the tendencies from integrating the model by one time step; 
- forcing option to add the prescribed forcing to the model  

see the changes in
- default src code: [spectral_dynamics.F90](../../../src/atmos_spectral/model/spectral_dynamics.F90)
- modified src code: [spectral_dynamics.F90](./src/spectral_dynamics.F90)  
diff ../../../src/atmos_spectral/model/spectral_dynamics.F90 ./src/spectral_dynamics.F90

Modified namelist parameters
- set the model truncation to zonally symmetric (i.e., fourier_inc = num_fourier)
- create climatology.nc using matlab; read the file as the initial condition (i.e., initial_state_option    = 'input'); set up the option in ic_from_external_file_nml for the input file
- adjust forcing and tendency values

See the changes in
- default namelist: [namelists](../default/run_ctl/namelists)
- namelist for run_ctl_forcing: [namelists](./run_ctl_forcing/namelists)
- namelist for run_ctl: [namelists](./run_ctl/namelists)  
diff ../default/run_ctl/namelists ./run_ctl/namelists

## Summary of directory contents
run_ctl_forcing/  
- model configuration to obtain the eddy forcing

run_ctl/  
- model configuration with prescribed eddy forcing

src/  
- pathnames of the source code 
