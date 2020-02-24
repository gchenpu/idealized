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

#------------------------
# control experiment
#------------------------
#model = fmspy.mo.FmsModel(name='hs_with_clouds', version="v3.0", exp_name="run_ctl_ts200_dh60", dayslist='( 6570 )', dt_atmos=1800, \
#                          queue_cmd="qsub", nnodes=1, ncores=16, nthreads=2)

#------------------------
# cooling the stratosphere
#------------------------
#model = fmspy.mo.FmsModel(name='hs_with_clouds', version="v3.0", exp_name="run_ctl_ts190_dh60", exp_name_src="run_ctl_ts200_dh60", dayslist='( 6570 )', dt_atmos=1800, \
#                          queue_cmd="qsub", nnodes=1, ncores=16, nthreads=2, \
#			  namelist_patch = {'hs_forcing_nml': {'t_strat': 190} })

#------------------------
# warming the troposphere
#------------------------
model = fmspy.mo.FmsModel(name='hs_with_clouds', version="v3.0", exp_name="run_ctl_ts200_dh60_3K", exp_name_src="run_ctl_ts200_dh60", dayslist='( 6570 )', dt_atmos=1800, \
                          queue_cmd="qsub", nnodes=1, ncores=16, nthreads=2, \
			  namelist_patch = {'hs_forcing_nml': {'t_zero': 318} })

#------------------------	  
# decrease the meridional temperature gradient
#------------------------
#model = fmspy.mo.FmsModel(name='hs_with_clouds', version="v3.0", exp_name="run_ctl_ts200_dh50", exp_name_src="run_ctl_ts200_dh60", dayslist='( 6570 )', dt_atmos=1800, \
#                          queue_cmd="qsub", nnodes=1, ncores=16, nthreads=2, \
#			  namelist_patch = {'hs_forcing_nml': {'delh': 50} })

#------------------------
# tropical upper tropospheric warming
#------------------------
#model = fmspy.mo.FmsModel(name='hs_with_clouds', version="v3.0", exp_name="run_ctl_ts200_dh60_TUW_A8_ph30", exp_name_src="run_ctl_ts200_dh60", dayslist='( 6570 )', dt_atmos=1800, \
#                          queue_cmd="qsub", nnodes=1, ncores=16, nthreads=2, \
#			  namelist_patch = {'hs_forcing_nml': {'newtonian': 'TUW',
#			                                        'AA': 8,
#								'ph0': 30} })

#model.setup2()        # run script is not executed
model.setup()        # run script is executed
