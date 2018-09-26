# Tutorial for running the FMS models
This is a tutorial for running a modified version of GFDL FMS models. This version has an option to override the advective zonal mean zonal wind (see <reference> for more details) in the dynamical core. The tutorial introduces a python script that compiles and run models automatically. 

The compilation and running of model is configured by a number of configuration files. Across different versions of model and run, only a small part of these config files differ. Therefore, a new configuration can be created by making edits on a batch of template config files. This is what this python script is mainly doing -- writing model/run config files by replacing content in templates according to the specifications. In addition to configure, the python script can also execute compilations and runs.
## A quick start for python script
A Python script should be executed at the root directory of the model. The script accesses a class called **FmsModel** defined in ```/fmspy/mo.py```. This class contains all the information (in the form of class variables) and available operations (in the form of class methods) for a specific model/run. The configuration is specified by passing parameters in the creation of an object of the class. The interface of the initialization function looks like this:
```python
def __init__(self, name="fullrad", version="default", exp_name="run_ctl", exp_name_src=[], platform="ncar", compile_opt="", echo_opt="unset echo", monthslist="(0)", dayslist="(1)", dt_atmos=600,  queue_cmd="", proj_code="UCOR0011", walltime="06:00:00", queue_name= "share", nnodes=1, ncores=1, nthreads=1, namelist_patch={}, grid_spec="", init_cond=""):
```
Although having many parameters, as will be shown later, only a small subset of the parameters need to be specified when doing certain things. Below is an example of how to use this class to compile and run the dynamical core with Held-Suarez forcing.
### Compile 
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
### Configure 
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
### Run 
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
### Debug
If a run crashes, the log can be found in ```/workdir/hs/v3.4e/delh70_ovrd_cntr/``` for the above example. The 4 configurations files (written by the Python script) that are actually read by the model can also be found in the same folder. 

