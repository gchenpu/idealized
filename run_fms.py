#===================================================================================
# main program to compile and run an experiment of one fms model
#
# need to specify user-defined output root_dir (default=../data) and input root_dir (default=../work) in fmspy/mo.py
# compile and run scripts are created in workdir (default=./workdir), which are temporary
# suggest to create symbolic links to these directories without modifying the python code
#===================================================================================
import fmspy.mo

# default values
# name="fullrad", version="default", exp_name="run_ctl", exp_name_src=[],                  
# platform="ncar", compile_opt="", echo_opt="unset echo", (monthslist, dayslist, secondslist) = ("(0)", "(1)", "(0)"), dt_atmos=600,            
# queue_cmd="", proj_code="UCOR0011", walltime="06:00:00", queue_name= "share", nnodes=1, ncores=1, nthreads=1
# namelist_patch={}, grid_spec="", init_cond="", model2plevel=False
#
# example of namelist_patch:
# namelist_patch = {'idealized_moist_phys_nml': {'radiation_scheme': 'two_stream'}, 'astronomy_nml': {'obliq': 23.5}}

#model = fmspy.mo.FmsModel(name="fullrad", dayslist='( 365 365 365 )', dt_atmos=600, \
#                          queue_cmd="qsub", walltime="06:00:00", nnodes=1, ncores=16, nthreads=2)

#model = fmspy.mo.FmsModel(name="gray",  dayslist='( 365 365 365 365 )', dt_atmos=600, \
#                          queue_cmd="qsub", nnodes=1, ncores=16, nthreads=2)

#model = fmspy.mo.FmsModel(name='hs',    dayslist='( 365 365 365 365 )', dt_atmos=1800, \
#                          queue_cmd="qsub", nnodes=1, ncores=16, nthreads=2)

model = fmspy.mo.FmsModel(name='martineau18', dayslist='( 120 )', dt_atmos=400, 
                          queue_cmd="qsub", nnodes=1, ncores=16, nthreads=2, model2plevel=True)

#model = fmspy.mo.FmsModel(name='hs_with_clouds', dayslist='( 1460 )', dt_atmos=1800, \
#                          queue_cmd="qsub", nnodes=1, ncores=16, nthreads=2)

#model = fmspy.mo.FmsModel(name='f90', version="spectral_tutorial", exp_name="run_ctl")

#model = fmspy.mo.FmsModel(name='f90', version="Tpdf", exp_name="run_ctl", queue_cmd="")

#model.setup2()        # run script is not executed
model.setup()        # run script is executed
