#===================================================================================
# main program to compile and run an experiment of one fms model
#
# need to specify user-defined INPUT_DIR (default=../work) and OUTPUT_DIR (default=../data) in fmspy/mo.py
# permanent compile and run scripts are available in MY_DIR (default=exp)
# temporatory compile and run scripts are saved in WORK_DIR (default=workdir)
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

#model = fmspy.mo.FmsModel(name='hs', version='wave_model', exp_name="run_ctl", dayslist='( 365 365 365 365 )', dt_atmos=1800, \
#                          queue_cmd="qsub", proj_code="UCLA0027", nnodes=1, ncores=16, nthreads=2)
			  
#model = fmspy.mo.FmsModel(name='hs', version='zonal_symmetric_model', exp_name="run_ctl_forcing", secondslist='( 1800 )', dt_atmos=1800, \
#                          queue_cmd="", proj_code="UCLA0027", nnodes=1, ncores=1, nthreads=1, \
#                          init_cond="/glade/u/home/chen/idealized/exp/hs/zonal_symmetric_model/run_ctl_forcing/climatology.nc")

#model = fmspy.mo.FmsModel(name='hs', version='zonal_symmetric_model', exp_name="run_ctl", dayslist='( 5 )', dt_atmos=1800, \
#                          queue_cmd="", proj_code="UCLA0027", nnodes=1, ncores=1, nthreads=1, \
#                          init_cond="/glade/u/home/chen/idealized/exp/hs/zonal_symmetric_model/run_ctl_forcing/climatology.nc")

model = fmspy.mo.FmsModel(name='hs', version='zonal_symmetric_model', exp_name="run_ctl", dayslist='( 365 365 )', dt_atmos=1800, \
                          queue_cmd="", proj_code="UCLA0027", nnodes=1, ncores=1, nthreads=1)

#model.setup2()        # run script is not executed
model.setup()        # run script is executed
