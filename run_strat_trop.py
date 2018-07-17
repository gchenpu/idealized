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

# control experiment
#model = fmspy.mo.FmsModel(name='strat_trop', exp_name="run_ctl_h0_ls-20_us20", dayslist='( 365 365 365 365 )', dt_atmos=600, 
#                          queue_cmd="qsub", nnodes=1, ncores=16, nthreads=1, model2plevel=True)

# create perturbation experiments using patch
# varying the vortex strength
#model = fmspy.mo.FmsModel(name='strat_trop', exp_name="run_ctl_h0_ls-30_us20", exp_name_src="run_ctl_h0_ls-20_us20", dayslist='( 365 365 365 365 )', dt_atmos=600, 
#                          queue_cmd="qsub", nnodes=1, ncores=16, nthreads=1, model2plevel=True, \
#                          namelist_patch = {'hs_forcing_nml': {'del_ts_nh': -30.}})

#model = fmspy.mo.FmsModel(name='strat_trop', exp_name="run_ctl_h0_ls-10_us20", exp_name_src="run_ctl_h0_ls-20_us20", dayslist='( 365 365 365 365 )', dt_atmos=600, 
#                          queue_cmd="qsub", nnodes=1, ncores=16, nthreads=1, model2plevel=True, \
#                          namelist_patch = {'hs_forcing_nml': {'del_ts_nh': -10.}})

# varying the topography height
#model = fmspy.mo.FmsModel(name='strat_trop', exp_name="run_ctl_h2w2_ls-20_us20", exp_name_src="run_ctl_h0_ls-20_us20", dayslist='( 365 365 365 365 )', dt_atmos=600, 
#                          queue_cmd="qsub", nnodes=1, ncores=16, nthreads=1, model2plevel=True, \
#                          namelist_patch = {'gaussian_topog_nml': {'topo_height': 2000., 
#			                                           'topo_wnum': 2} })

#model = fmspy.mo.FmsModel(name='strat_trop', exp_name="run_ctl_h4w2_ls-20_us20", exp_name_src="run_ctl_h0_ls-20_us20", dayslist='( 365 365 365 365 )', dt_atmos=600, 
#                          queue_cmd="qsub", nnodes=1, ncores=16, nthreads=1, model2plevel=True, \
#                          namelist_patch = {'gaussian_topog_nml': {'topo_height': 4000., 
# 			                                            'topo_wnum': 2} })

model = fmspy.mo.FmsModel(name='strat_trop', exp_name="run_ctl_h6w2_ls-20_us20", exp_name_src="run_ctl_h0_ls-20_us20", dayslist='( 365 365 365 365 )', dt_atmos=600, 
                          queue_cmd="qsub", nnodes=1, ncores=16, nthreads=1, model2plevel=True, \
                          namelist_patch = {'gaussian_topog_nml': {'topo_height': 6000., 
 			                                            'topo_wnum': 2} })

#model.setup2()        # run script is not executed
model.setup()        # run script is executed
