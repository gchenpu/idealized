# Tutorial for running the FMS models
This is a tutorial for running a modified version of GFDL FMS models. This version has an option to override the advective zonal mean zonal wind (see <reference> for more details) in the dynamical core. The first part of the tutorial introduces a python script that compiles and run models automatically. The second part introduces the modifications to the model code. The third part document the setups for a few standrad and overriding simulations.


## Part I: A quick start (introducing the python script)
The compilation and running of model is configured by a number of configuration files. Across different versions of model and run, only a small part of these config files differ. Therefore, a new configuration can be created by making edits on a batch of template config files. This is what this python script is mainly doing -- writing model/run config files by replacing content in templates according to the specifications. In addition to configure, the python script can also execute compilations and runs.

The Python script should be executed at the root directory of the model. The script accesses a class called **FmsModel** defined in ```/fmspy/mo.py```. This class contains all the information (in the form of class variables) and available operations (in the form of class methods) for a specific model/run. The configuration is specified by passing parameters in the creation of an object of the class. The interface of the initialization function looks like this:
```python
def __init__(self, name="fullrad", version="default", exp_name="run_ctl", exp_name_src=[], platform="ncar", compile_opt="", echo_opt="unset echo", monthslist="(0)", dayslist="(1)", dt_atmos=600,  queue_cmd="", proj_code="UCOR0011", walltime="06:00:00", queue_name= "share", nnodes=1, ncores=1, nthreads=1, namelist_patch={}, grid_spec="", init_cond=""):
```
Although having many parameters, as will be shown later, only a small subset of the parameters need to be specified when doing certain things. Below is an example of how to use this class to compile and run the dynamical core with Held-Suarez forcing.
### 1. Compile 
A minimalist script that compile a version of the dynamical core would look like this:
```python
import fmspy.mo
model = fmspy.mo.FmsModel(name='hs', version='v3.4e')
model.setup_mkmf()
model.setup_compile()
model.exec_compile()
```
Here, an object called _model_ is an instance of the _FmsModel_ class. It is initialized with two parameters specified: _name_ and _version_. The _name_ is set to _"hs"_ which corresponds to Held-Suarez setup. One can choose other options like _"gray"_ or _"full_rad"_. This tells the object where to find the source code for a specific model. 

The _version_ is set to a name defined by the user. This tells the model where to find the modified source code and include them in the compilation. Here it is set to _"v3.4e"_. This means in the folder ```/exp/hs/v3.4e/src/```, there should be modified source code for this specific version. A peek into this folder would look like this:
![logo](https://raw.githubusercontent.com/dkwovn/fmspy_tutorial/master/v34e_src.png?token=AfyR3eomXh5aqgb1jD51RKPD4O_GsC2Xks5bkwsgwA%3D%3D)
Note that there is a file called ```path_name```. This should be copied from the ```/exp/hs/default/src``` folder. Also note that original versions of these F90 files exist in the ```/src/``` directory and will be ignored in the compilation. 

If you need to change the source code of the model and make a new version, you can create a folder with the new version name ```/exp/hs/new_version_name/src```, copy the ```path_name``` into this folder and the F90 files you are willing to modify. Then modify those F90 files. Lastly, execute a compile script like the above one except set
```python
model = fmspy.mo.FmsModel(name = 'hs', version = 'new_version_name')
```
in the second line.
### 2. Configure 
With a compiled model, we can run experiments using various configurations. Without the _FmsModel_ class, this would usually involve editing the following 4 files.
1. namelists
   A template of default namelists is stored as ```/exp/hs/default/run_ctl/namelists```. It is usually the most frequently changed file as it specifies almost all the running configurations such as resolution, damping, physics etc. The overriding options are also configured through namelists.
2. diag_table
   A template of default diag_table is stored as ```/exp/hs/default/run_ctl/diag_table```. It specifies the variables to output, as well as the frequency and nature of the output (instantaneous or average). 
3. field_table
    A template of field_table is stored as ```/exp/hs/default/run_ctl/field_table```. It specifies all the information related to individual tracers in the model.
4. runscript
   A template of runscript is stored as ```/bin/run.template```. This template is not a working runscript but a meta data for _FmsModel_ class to fill in the parameters on the fly. 

The _FmsModel_ object provides methods to modify/create only (1) namelists and (4) runscript. The (2) diag_table and (3) field_table used in the model can be specified in the runscript but cannot be changed on the fly as (1) and (4).

A typical configuration of _FmsModel_ looks like this:
```Python
import fmspy.mo
monthslist = [0 0 0 0]
dayslist = [365 365 365 365]
template_path = '/glade/u/home/lchengji/fms/idealized/exp/hs/default/run_ctl'
name_dict = {}
name_dict['hs_forcing_nml'] = {'delh': 70}
name_dict['spectral_dynamics_nml'] = \
   {'ovrd_option': 'full',\
    'ovrd_folder': '/glade/u/home/lchengji/data/hs/runs_v3.3/control1',\
    'ovrd_file_prefix': 'hs.6hourly.',\
    'ovrd_file_suffix': '.nc',\
    'ovrd_freq_option': '6hourly'}

model = fmspy.mo.FmsModel(name='hs', version='v3.4e', exp_name='delh70_ovrd_cntr', exp_name_src=template_path, monthslist=monthslist, dayslist=dayslist, namelist_patch=name_dict)
```
All the variable assignments before the instantiation is just to make the code a bit cleaner. Let's first focus on the actual instantiation at the last line. The parameters _name_ and _version_ tell model what executable to use. So make sure this version of model has been compiled according to the steps in the previous section. The _exp_name_ specifies the name of run folder that will be created in ```/exp/hs/v3.4e/```. _exp_name_src_ tells model where to find the templates for namelists, diag_table, and field_table. The rest of the parameters either configure the namelists or the runscript.
##### Namelists 
All the modifications of the namelists is wrapped into a single parameter called _namelist_patch_. It is supposed to be a Python dictionary, i.e. nested key-value pairs. The dictionary for this run is constructed in the previous lines and is called _name_dict_. In this example, the _name_dict_ has two keys: 'hs_forcing_nml' and 'spectral_dynamics_nml', which are the names of two model nameslists. The values corresponding to these two keys are two additionaly dictionaries specifying the two namelists respectively. All the namelists values are case insensitive. 

All the variables in the 'spectral_dynamics_nml' are unique to this version of model and configure the behavior of the overriding. The 'ovrd_option' can be 'full' meaning it's fully on, or 'none' meaning it's off. In the case of 'none', all the following variables would be ignored. The 'ovrd_folder' tells model where to look for the data storing the overriding zonal wind. 'ovrd_file_prefix' and 'ovrd_file_suffix' tells model the form of these overriding file names. It is assumed to be "prefix + year + suffix". The last variable 'ovrd_freq_option' specifies the frequency at which the overriding wind is updated. 
##### Runscript
The two parameters that involve runscript changes are monthslist (might have been obsolete) and dayslist. The dayslist is a Python list that listed the number of days at each year. In this example, the run will cover 4 years with 365 days in each year. 
### 3. Run 
Before we want to hit the go button, there are a few more options we might want to configure. These options don't affect the model output but affect performance, cost etc.
```Python
model = fmspy.mo.FmsModel(..., dt_atmos=1800, queue_name="share", queue_cmd="qsub", nnodes=1, ncores=16, nthreads=2)
```
These options are implemented by editing the runscript. _dt_atmos_ configure the time step length for integration. Make sure it satisfies the CFL criteria, especially when you are changing spatial resolution. The _queue_name_ specifies the queue you will submit to. The cost for different queues can be found here (https://bit.ly/2MRl2kK). _queue_cmd_ specifies the command you will use to submit the job. _nnodes_ and _ncores_ specify the number of nodes and cores used in the job. It is recommanded to refer to the link above as every queue has a limit for how many cores can be used and how the charge is calculated etc.

Once all these things are set up, you can use two methods to write the configuration files and excute the run respectively:
```Python
model.setup_run()
model.exec_run()
```
To sum up, a Python script to configure and run a job would look like:
```Python
import fmspy.mo

monthslist = [0 0 0 0]
dayslis = [365 365 365 365]
template_path = '/glade/u/home/lchengji/fms/idealized/exp/hs/default/run_ctl'
name_dict = {}
name_dict['hs_forcing_nml'] = {'delh': 70}
name_dict['spectral_dynamics_nml'] = \
   {'ovrd_option': 'full',\
    'ovrd_folder': '/glade/u/home/lchengji/data/hs/runs_v3.3/control1',\
    'ovrd_file_prefix': 'hs.6hourly.',\
    'ovrd_file_suffix': '.nc',\
    'ovrd_freq_option': '6hourly'}

model = fmspy.mo.FmsModel(name='hs', version='v3.4e', exp_name='delh70_ovrd_cntr', exp_name_src=template_path, monthslist=monthslist, dayslist=dayslist, namelist_patch=name_dict, queue_name="share", queue_cmd="qsub", nnodes=1, ncores=16, nthreads=2)

model.setup_run()
model.exec_run()
```
### 4. Debug
If a run crashes, the log can be found in ```/workdir/hs/v3.4e/delh70_ovrd_cntr/``` for the above example. The 4 configurations files (written by the Python script) that are actually read by the model can also be found in the same folder. 

## Part II: Modifications for FMS dynamical core
### spectral_dynamics.F90
__The version for illustration is Held-Suarez model ```v3.4f```. The source code is in ```/glade/u/home/lchengji/fms/idealized/exp/hs/v3.4f/src```. Most of the modifications have comments starting with ```cliu30/cliu31/cliu33```.__ The goal is to calculate the differential advection term on the RHS of (17) in the manscript.
#### 1. Namelist variables
The overriding in the dynamical core is optional. The switch that turns on/off the overriding is a namelist variable:
```Fortran
character(len=200) :: ovrd_option      = 'none' ! The overriding option, can be 'full' or 'ubt' (barotropic)
```
When ```ovrd_option``` is set to be ```'none'```, the overriding is turned off and the model is standard. When it is set to be ```'ubt'```, only the barotropic component of zonal-mean zonal wind will be overrode. When it is set to be ```'full'```, the full zonal-mean zonal wind would be overrode.

For the ```'full'``` overriding, the overriding wind consists of two parts. The first part is a perturbation field that doesn't change over time. This part is optional and is configured in the ```spectral_dynamics_init(...)``` subroutine. The second part of the overriding is a 6-hourly field read from another simulation. This part is configured in the ```spectral_dynamics(...)``` subroutine.
#### 1st part (invariant perturbation)
The type of perturbation is set by a namelist variable:
```Fortran
character(len=200) :: ovrd_pert_option = 'none' ! whether to add a perturbation to the overriding field or not. Can be 'const', 'gauss', or 'input'
```
When ```ovrd_pert_option = 'gauss'```, the parameters for the Gaussian perturbation are set by these namelist variables:
```Fortran
! parameters for ovrd_pert_option = 'gauss' / 'const'	
real :: lato_center         = 0.0   ! the center latitude for Gaussian perturbation
real :: lato_width          = 10.0  ! the meridional width of the Gaussian perturbation
real :: sigo_center         = 0.5   ! the center vertical level
real :: sigo_width          = 0.075 ! the vertical width
real :: uo_pert             = 4.0   ! the amplitude of the Gaussian perturbation
                                    ! also used as the amplitude of the perturbation for ovrd_pert_option = 'const'
```
When ```ovrd_pert_option = 'const'```, the constant perturbation is set by:
```Fortran
real :: uo_pert             = 4.0   ! the amplitude of the Gaussian perturbation
                                    ! also used as the amplitude of the perturbation for ovrd_pert_option = 'const'
```
When ```ovrd_pert_option = 'input'```, it means the perturbation is read from the first time step of a netCDF file. The file name is set by:
```Fortran
character(len=200) :: ovrd_pfile_name  = '' ! the perturbation file name for ovrd_pert_option = 'input'
```
##### 2nd part (6hourly input)
For the part that reads in data 6-hourly, the file information is set by these variables:
```Fortran
character(len=200) :: ovrd_freq_option = '6hourly' ! the frequency of overriding, currently only support 6hourly
character(len=200) :: ovrd_folder      = '' ! the folder path for (6-hourly) overriding files
character(len=200) :: ovrd_file_prefix = 'fullrad.' ! the prefix for overriding files
character(len=200) :: ovrd_file_suffix = '.atmos_4xday.nc' ! the suffix of overriding files
character(len=200) :: ovrd_file_name   = '' ! this is generated by the model as prefix + yr_str + suffix
```
#### 2. Work variables
The following variables are declared to do the actual work of overriding. Their meanings are indicated in the comments. Note that some obsolete variables are still retained in the source code but not listed here.
```Fortran
real, allocatable, dimension(:, :, :) :: ug_diff  ! the gridded difference between specified overriding u wind and the current u wind produced by the model
real, allocatable, dimension(:, :, :) :: ug_pert  ! the gridded perturbation field

complex, allocatable, dimension(:, :, :) :: vors_diff   ! the spectral difference between the vorticity of specified u wind and the current vorticity produced by the model
complex, allocatable, dimension(:, :, :) :: dt_vors_tmp ! a temperary variable to store vors tendency

real, allocatable, dimension(:, :, :) :: u_spcf      ! gridded field of current specified u wind
real, allocatable, dimension(:, :, :) :: u_spcf_prev ! gridded field of previous 6-hourly u wind
real, allocatable, dimension(:, :, :) :: us_zm       ! zonal mean of specified u wind
real, allocatable, dimension(:, :, :) :: dt_vorg_tmp ! a temperary variable to store vorg
integer :: id_us_zm ! id for outputting us_zm field
```
#### 3. Configure perturbation part
The code that setup the perturbation is in:
```Fortran
subroutine spectral_dynamics_init(Time, Time_step_in, tracer_attributes, dry_model_out, nhum_out, ocean_mask, n_days, n_years)
```
The block of code reads:
```Fortran
!-- configure the perturbation (not a function of time) --
!- Gaussian -
if (uppercase(trim(ovrd_pert_option)) == 'GAUSS') then
	allocate(temp_full(num_levels))
	allocate(sig1d(num_levels))
	allocate(temp_half(num_levels + 1))
    call get_deg_lat(dg_lat)
	call pressure_variables(temp_half, temp_half, sig1d, temp_full, reference_sea_level_press)
	u_spcf = zonal_mean(ug(:, :, :, 1))
	sig1d = sig1d / reference_sea_level_press
  	do k = 1, num_levels
  	    do j = js, je
   		    ug_pert(:, j, k) = uo_pert * exp(-((dg_lat(j) - lato_center) / lato_width)**2 &
    			                   -((sig1d(k) - sigo_center) / sigo_width)**2)
   		    ug_pert(:, j, k) = ug_pert(:, j, k) + uo_pert * exp(-((dg_lat(j) + lato_center) / lato_width)**2 &
    			                   -((sig1d(k) - sigo_center) / sigo_width)**2)
   	    end do
   	end do
	deallocate(temp_full)
	deallocate(temp_half)
	deallocate(sig1d)

!- Specified by an input file -
else if (uppercase(trim(ovrd_pert_option)) == 'INPUT') then
	call read_data(trim(ovrd_pfile_name), 'diff', ug_pert(:, :, :), grid_domain, timelevel = 1)
	!note that the variable name is assumed to be 'diff'
	ug_pert = zonal_mean(ug_pert)

!- Constant -
else if (uppercase(trim(ovrd_pert_option)) == 'CONST') then
	ug_pert(:, :, :) = uo_pert

!- No perturbation -
else if (uppercase(trim(ovrd_pert_option)) == 'NONE') then
	ug_pert = 0
else 
	call error_mesg('spectral_dynamics', trim(ovrd_pert_option)//'is not valid!', FATAL)
end if
```
Note that for ```ovrd_pert_option = 'input'```, the variable name in the input file must be ```diff```.
#### 4. Configure the 2nd part
The input file for overriding zonal-mean zonal wind contains 6-hourly wind fields. However, the model time step is much smaller than 6 hour. So at each model time step, the overriding wind is obtained by linear interpolation using two consecutive 6-hourly input ```u_spcf_prev``` and ```u_spcf```. 
##### Initialization of ```u_spcf```
In the:
```Fortran
subroutine spectral_dynamics_init(Time, Time_step_in, tracer_attributes, dry_model_out, nhum_out, ocean_mask, n_days, n_years)
```
The ```u_spcf``` is initialized as:
```Fortran
!-- initialize u_spcf --
if (uppercase(trim(ovrd_option)) /= 'NONE') then
	write(year_str, "(i5)") curr_year
	!-- construct ovrd_file_name as prefix + year + suffix --
    ovrd_file_name = trim(ovrd_folder) // '/' // trim(ovrd_file_prefix) // & 
		trim(adjustl(year_str)) // trim(ovrd_file_suffix)

	if (file_exist(trim(ovrd_file_name))) then
	    call read_data(trim(ovrd_file_name), 'ucomp', u_spcf(:, :, :), grid_domain, timelevel = 1)
    else
		call error_mesg('spectral_dynamics', trim(ovrd_file_name)//' does not exist!', FATAL)
	end if
end if
```
##### Reading ```u_spcf```
After the initialization, every 6 hour, the old ```u_spcf``` is assigned to ```u_spcf_prev``` and the new ```u_spcf``` is read in 
```Fortran
subroutine spectral_dynamics(Time, psg_final, ug_final, vg_final, tg_final, tracer_attributes, grid_tracers_final, &
                             time_level_out, dt_psg, dt_ug, dt_vg, dt_tg, dt_tracers, wg_full, p_full, p_half, z_full)
```
The code reads:
```Fortran
call get_time(Time, seconds, days)
	if (uppercase(trim(ovrd_freq_option)) == '6HOURLY') then
		if(mod(seconds, time_int) == 0 .and. file_exist(trim(ovrd_file_name))) then
		! time_int = 6 hour (6 * 60 * 60 seconds)
			u_spcf_prev = u_spcf
			time_step = mod(max(days, 0), num_days) * 4 + seconds/time_int + 1
			call read_data(trim(ovrd_file_name), 'ucomp', u_spcf(:, :, :), grid_domain, timelevel = time_step)
		else if (.not.file_exist(trim(ovrd_file_name))) then
			call error_mesg('spectral_dynamics', trim(ovrd_file_name)//' does not exist!', FATAL)
		end if
	else 
	    ...
	end if
```
Other options of reading overriding field can be added in the ```else``` block. 
#### 5. Calculating ```ug_diff```
```ovrd_option``` determines the way to calculate ```ug_diff```. For ```ovrd_option = 'ubt'``` (barotropic), the code reads:
```Fortran
    else if(uppercase(trim(ovrd_option)) == 'UBT') then
		us_zm = zonal_mean_vm(u_spcf)
		ug_diff = zonal_mean_vm(time_interp(u_spcf_prev, u_spcf, seconds) - ug(:, :, :, current))
		vors_diff = 0
```
For ```ovrd_option = 'full'```, the code reads:
```Fortran
    else if (uppercase(trim(ovrd_option)) == 'FULL') then
		us_zm = zonal_mean(u_spcf) + ug_pert
		ug_diff = zonal_mean(time_interp(u_spcf_prev, u_spcf, seconds) - ug(:, :, :, current)) + ug_pert
		vors_diff = 0
```
#### 6. Adding the differential advection term
After ```ug_diff``` is calculated, the differential advection term is added to the tendency of temperature and vorticity respectively (see equation (17) in the manuscript):
```Fortran
!--- temperature advection ---
call horizontal_advection(ts(:,:,:,current), ug_diff, vg(:,:,:,current) * 0, dt_tg_tmp)
dt_tg_tmp = dt_tg_tmp + tg_diff
!==============
call trans_grid_to_spherical(dt_tg_tmp, dt_ts)
```
```Fortran
dt_vorg_tmp = 0.
!--- vorticity advection ---
call horizontal_advection(vors(:, :, :, current), ug_diff, ug_diff * 0, dt_vorg_tmp)
call trans_grid_to_spherical(dt_vorg_tmp, dt_vors_tmp)
dt_vors = dt_vors + dt_vors_tmp
```
This concludes the modifications for overriding. 
## Part III: Climate change runs
Below are a few examples for setting up climate change simulations in three models of different hierarchy. 
### 1. Held-Suarez model
There a few ways to change the thermal or mechanical forcing in the model. To change the horizontal gradient of the equilibruim temperature, one can alter a namelist variable ```delh``` that specify the equator-to-pole temperature difference. To set this parameter to be 70K, the python script would look like this:
```Python
import fmspy.mo
...
name_dict = {}
name_dict['hs_forcing_nml'] = {'delh': 70}
model = fmspy.mo.FmsModel(name='hs', version='v3.4e', exp_name='delh70_ovrd_cntr', exp_name_src=template_path, monthslist=monthslist, dayslist=dayslist, namelist_patch=name_dict)
```
Another parameter is stratosphere temperature. To change it, just replace the fourth row of the above example with:
```Python
name_dict['hs_forcing_nml'] = {'t_strat': 196}
```
This sets the stratosphere temperature to be 196K. 

To change the friction strength for the zonal-mean zonal wind in the boundary layer, the namelist variables are:
```Python
name_dict['hs_forcing_nml'] = {'rayleigh': 'eddy', 'kf_z': -1.5}
```
The damping e-folding time scale has been set to 1.5 days. 
### 2. Gray radiation model
To mimic global warming/cooling in the Gray radiation model, one can change the opitical depth for longwave radiation. This can be set via variables in the namelist ```two_stream_gray_rad_nml``` in a similar way in the python script:
```Python
name_dict['two_stream_gray_rad_nml'] = {'ir_tau_pole': 1.8, 'ir_tau_eq': 7.2}
```
The above code set the optical depth to be 1.8 at the pole and 7.2 at the equator. 

Note that the integration time step is decreased to 600s in these simulations to avoid numerical instability:
```Python
 model = fmspy.mo.FmsModel(name='gray', version='v3.4d', ..., dt_atmos=600, ...)
```
### 3. Full radiation model
The warming in the full radiation model is simulated by increase the constant CO2 concentration. For example:
```Python
name_dict['full_radiation_driver_nml'] = {'sw_rad_time_step': 1200, 'rco2': 0.00074}
```
In the above example, the CO2 concentration is set to be 740 ppm, about twice as in the control simulation. Also, the integration time step for shortwave radiation is increased to 1200s to speed up the run.

The integration time step for the dynamical core is 600s, same as the gray radiation model.
```Python
 model = fmspy.mo.FmsModel(name='fullrad', version='v3.3', ..., dt_atmos=600, ...)
```
