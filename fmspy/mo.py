#===================================================================================
# mo.py: module to handle mkmf, compile, and run scripts in the fmspy package
# default user-defined dir of files:  "./exp"
# default user-defined output root_dir: "../data"
# default user-defined input root_dir: "../work"
#===================================================================================
import os
import sys
import shutil
import subprocess
import tarfile
import fmspy.io

class FmsModel:

    #===================================================================================
    # initialize class variables
    #===================================================================================
    def __init__(self, name="fullrad", version="default", exp_name="run_ctl", exp_name_src=[],                        \
        platform="ncar", compile_opt="", echo_opt="unset echo", monthslist="(0)", dayslist="(1)", dt_atmos=600,       \
	queue_cmd="", proj_code="UCOR0011", walltime="06:00:00", queue_name= "share", nnodes=1, ncores=1, nthreads=1, \
	namelist_patch={}, grid_spec="", init_cond=""):
        # example of namelist_patch when exp_name_src is not []
        # namelist_patch = {'idealized_moist_phys_nml': {'radiation_scheme': 'two_stream'}, 'astronomy_nml': {'obliq': 23.5}}

        if not name in ["fullrad", "gray", "hs", "hs_with_clouds", "f90"]:
            raise Exception("\"%s\" is not a model name!" % name)

        self.name            = name                                                 # model name
        self.version         = version                                              # code version
        self.exp_name        = exp_name                                             # exp name
        self.exp_name_src    = exp_name_src                                         # src exp used to create new exp

        self.platform        = platform                                             # platform
        self.compile_opt     = compile_opt                                          # compile option: "" for regular run and "debug" for debug
        self.echo_opt        = echo_opt                                             # unset echo or set echo in compile and run shell scripts        

        self.monthslist      = monthslist
        self.dayslist        = dayslist
        self.dt_atmos        = dt_atmos                                             # model time step in seconds

        self.queue_cmd       = queue_cmd                                            # "" or "qsub"
        self.proj_code       = proj_code
        self.walltime        = walltime
        self.queue_name      = queue_name
        self.nnodes          = nnodes
        self.ncores          = ncores
        self.nthreads        = nthreads

        self.namelist_patch  = namelist_patch
        self.grid_spec       = grid_spec
        self.init_cond       = init_cond

        self.mydir           = "exp"                                                ### user-defined dir of files
        self.workdir         = "workdir"                                            ### user-defined work dir
        self.cwd             = os.getcwd()
        self.output_root     = os.path.abspath(os.path.join(self.cwd, "../data"))   ### user-defined output root_dir
        self.input_root      = os.path.abspath(os.path.join(self.cwd, "../work"))   ### user-defined input root_dir

        if (self.exp_name_src != []) and (self.exp_name_src != self.exp_name):
            src              = os.path.join(".", self.mydir, self.name, self.version, self.exp_name_src)
            dst              = os.path.join(".", self.mydir, self.name, self.version, self.exp_name)
            if(not os.path.exists(src)):
                raise Exception("%s does not exist!" % src)
            else:
                print("\nSETTING UP \"%s\"" % dst)
            fmspy.io.copy_exp_setup(src, dst, namelist_patch=self.namelist_patch, option='silent')

        # set the paths for mkmf template, compile, run scripts
        if self.compile_opt == "":
            self.mkmf_filename    = os.path.join(".", self.workdir, "bin", "mkmf.template."+self.platform)
            self.compile_filename = os.path.join(".", self.workdir, self.name, "compile."+self.version)
            self.run_filename     = os.path.join(".", self.workdir, self.name, self.version, self.exp_name, "run")
        elif self.compile_opt == "debug":
            self.mkmf_filename    = os.path.join(".", self.workdir, "bin", "mkmf.template."+self.platform+"."+self.compile_opt)
            self.compile_filename = os.path.join(".", self.workdir, self.name, "compile."+self.version+"."+self.compile_opt)
            self.run_filename     = os.path.join(".", self.workdir, self.name, self.version+"."+self.compile_opt, self.exp_name, "run")

        # platform dependent configurations
        if self.platform == "ncar":
            # compile options in mkmf template
            self.FC          = "mpif90"
            self.LD          = "mpif90"
            self.CC          = "mpicc"
            self.OPENMP      = "-qopenmp"
            if self.compile_opt == "":
                self.FFLAGS = "$(CPPFLAGS) %s -fno-alias -stack_temps -safe_cray_ptr -ftz -assume byterecl -g -O2 -i4 -r8 -nowarn -Wp,-w" % (self.OPENMP)
            elif self.compile_opt == "debug":
                self.FFLAGS = "$(CPPFLAGS) %s -fltconsistency -stack_temps -safe_cray_ptr -ftz -assume byterecl -g -O0 -i4 -r8 -check -check noarg_temp_created -warn noerrors -debug variable_locations -traceback" % (self.OPENMP)

            # compiler info in compile and run scripts
            self.compiler_info= \
"""## set up the intel compiler and mpi
source $MODULESHOME/init/csh
module purge
module load intel
module load impi"""

            # queue info in run script
            self.queue_info  = \
"""#PBS -N %(name)s
#PBS -A %(proj_code)s
#PBS -l walltime=%(walltime)s
#PBS -q %(queue_name)s
#PBS -l select=%(nnodes)s:ncpus=%(ncores)s:mpiprocs=%(ncores)s:ompthreads=%(nthreads)s
#PBS -j oe
#PBS -o %(output_dir)s/output
#PBS -m abe """  % {'name':self.name, 'proj_code':self.proj_code, 'walltime':self.walltime, 'queue_name':self.queue_name, \
                         'nnodes':self.nnodes, 'ncores':self.ncores, 'nthreads':self.nthreads,                                 \
                         'output_dir':os.path.dirname(os.path.abspath(self.run_filename))}

    #===================================================================================
    # setup the mkmf directory
    #===================================================================================
    def setup_mkmf(self):
        # create the mkmf directory
        mkmf_dir = os.path.dirname(self.mkmf_filename)
        print("\nSETTING UP \"%s\"" % mkmf_dir)
        if (not os.path.exists(mkmf_dir)):
            os.makedirs(mkmf_dir)

        # copy bin files
        src      = os.path.join(".", "bin")
        dst      = mkmf_dir
        fmspy.io.copyfolder(src, dst, files_to_exclude=["mkmf.template", "compile.template", "run.template"], option='silent')

        # create local mkmf template
        output_format = {'cwd':self.cwd, 'FC':self.FC, 'LD':self.LD, 'CC':self.CC, 'FFLAGS':self.FFLAGS, 'OPENMP':self.OPENMP}
        fmspy.io.write_script(self.mkmf_filename, os.path.join(".", "bin/mkmf.template"), output_format)

    #===================================================================================
    # setup the compile directory
    #===================================================================================
    def setup_compile(self):
        # create the compile directory
        compile_dir = os.path.dirname(self.compile_filename)
        print("\nSETTING UP \"%s\"" % compile_dir)
        if (not os.path.exists(compile_dir)):
            os.makedirs(compile_dir)

        # copy path_names
        full_src_name = os.path.join(".", self.mydir, self.name, self.version, "src/path_names")
        full_dst_name = os.path.join(compile_dir, "path_names")
        fmspy.io.copyfile(full_src_name, full_dst_name, option="overwrite")

        # copy my src code of the model to the src directory (remove existing dst folder before copying)
        src      = os.path.join(".", self.mydir, self.name, self.version, "src")
        src_files= fmspy.io.get_src_files(src, fmspy.io.SRC_TYPES)
        if src_files:
            dst      = os.path.join("./src", "my", self.name, self.version)
            print("\nSETTING UP \"%s\"" % dst)
            if (os.path.exists(dst)):
                shutil.rmtree(dst)
            fmspy.io.copyfolder2(src, src_files, dst, option="overwrite")
            myList   = '(\"' + " ".join('$sourcedir/%s/$name/$version/' % "my" +file for file in src_files) + '\")'
        else:
            myList  ='(\"\")'
        #print(myList)

        # create compile script
        if self.name=="gray" or self.name=="fullrad":
            moist_model_opt = "-DIDEALIZED_MOIST_MODEL"
        else:
            moist_model_opt = ""

        output_format = {'echo_opt':self.echo_opt, 'compile_opt':self.compile_opt, 'name':self.name, 'version':self.version, \
                         'platform':self.platform, 'myList':myList, 'compiler_info':self.compiler_info, "moist_model_opt":moist_model_opt}
        fmspy.io.write_script(self.compile_filename, os.path.join(".", "bin/compile.template"), output_format, option='overwrite')

    #===================================================================================
    # execute the compile script
    #===================================================================================
    def exec_compile(self):
        # exec the compile script
        print("\nCOMPILING \"%s\"" % self.compile_filename)
        os.chdir(os.path.dirname(self.compile_filename))

        #os.system('./' + os.path.basename(self.compile_filename))
        process = subprocess.Popen('./' + os.path.basename(self.compile_filename), stdout=subprocess.PIPE)
        std     = process.stdout.read().decode("utf-8")
        #print(std)
	
        os.chdir(self.cwd)

        if 'successful' in std.split("\n")[-2]:
            return 0
        else:
            return 1

    #===================================================================================
    # setup the input Data and return the pathname of inputData
    #===================================================================================
    def setup_inputData(self):
        src       = os.path.join(".", self.mydir, self.name, "input_data")
        if os.path.exists(src):
            src_files = fmspy.io.get_src_files(src, fmspy.io.DAT_TYPES)
        else:
            return '\"\"'

        if src_files:
            dst = os.path.join(self.input_root, self.name, "input_data")
            if (os.path.exists(dst)):
                print("\"%s\" already exists!" % dst)
            else:
                print("SETTING UP: \"%s\"" % dst)
                os.makedirs(dst)
                fmspy.io.copyfolder2(src, src_files, dst)
                download_inputData = True

            os.chdir(dst)
            if os.path.exists("./inputData.csh"):
                if download_inputData:
                    os.system("./inputData.csh")
                os.remove('inputData.csh')
            for fname in os.listdir("."):
                if (fname.endswith("tar.gz")):
                    tar = tarfile.open(fname, "r:gz")
                    tar.extractall()
                    tar.close()
                    os.remove(fname)
                elif (fname.endswith("tar")):
                    tar = tarfile.open(fname, "r:")
                    tar.extractall()
                    tar.close()
                    os.remove(fname)
            os.chdir(self.cwd)
            return "(%s/$name/input_data/*)" % self.input_root
        else:
            return '\"\"'

    #===================================================================================
    # setup the run directory
    #===================================================================================
    def setup_run(self):
        # create the run directory (removing existing run directory first)
        run_dir = os.path.dirname(self.run_filename)
        print("\nSETTING UP \"%s\"" % run_dir)
        if (os.path.exists(run_dir)):
            shutil.rmtree(run_dir)
        os.makedirs(run_dir)
		
        # copy diag_table, data_table, field_table, and namelists
        src      = os.path.join(".", self.mydir, self.name, self.version, self.exp_name)
        dst      = run_dir

        fmspy.io.copy_exp_setup(src, dst, option='overwrite')

        if os.path.exists(dst+'/diag_table'):
            diagTable        = "$rootdir/$name/$version/$exp_name/diag_table"
            diagTable2       = "$rootdir/$name/$version.$compile_opt/$exp_name/diag_table"
        else:
            diagTable        = "\"\""
            diagTable2       = "\"\"" 

        if os.path.exists(dst+'/data_table'):
            dataTable        = "$rootdir/$name/$version/$exp_name/data_table"
            dataTable2       = "$rootdir/$name/$version.$compile_opt/$exp_name/data_table"
        else:
            dataTable        = "\"\""
            dataTable2       = "\"\"" 

        if os.path.exists(dst+'/field_table'):
            fieldTable        = "$rootdir/$name/$version/$exp_name/field_table"
            fieldTable2       = "$rootdir/$name/$version.$compile_opt/$exp_name/field_table"
        else:
            fieldTable        = "\"\""
            fieldTable2       = "\"\"" 

        if os.path.exists(self.grid_spec):
            fmspy.io.copyfile(self.grid_spec, dst+"/"+os.path.basename(self.grid_spec), option='overwrite')
            gridSpec        = "$rootdir/$name/$version/$exp_name/"+os.path.basename(self.grid_spec)
            gridSpec2       = "$rootdir/$name/$version.$compile_opt/$exp_name/"+os.path.basename(self.grid_spec)
        else:
            gridSpec        = "\"\""
            gridSpec2       = "\"\"" 

        if os.path.exists(self.init_cond):
            fmspy.io.copyfile(self.init_cond, dst+"/"+os.path.basename(self.init_cond), option='overwrite')
            initCond        = "$rootdir/$name/$version/$exp_name/"+os.path.basename(self.init_cond)
            initCond2       = "$rootdir/$name/$version.$compile_opt/$exp_name/"+os.path.basename(self.init_cond)
        else:
            initCond        = "\"\""
            initCond2       = "\"\"" 

        # create run script
        if self.name == 'gray' or self.name == 'fullrad':
            calendar  = 'noleap'
        else:
            calendar  = 'no_calendar'

        inputData        = self.setup_inputData()
        FTMPDIR          = "%s/tmp/%s" % (self.output_root, "$HOST.`date '+%m-%d@%H:%M:%S'`.id$$")
        other_time_stamp = "tmp`date '+%j%H%M%S'`"

        output_format = {'queue_info':self.queue_info, 'echo_opt':self.echo_opt, 'compile_opt':self.compile_opt, 'cwd':self.cwd, 'workdir':self.workdir, \
                         'name':self.name, 'version':self.version, 'exp_name':self.exp_name, 'FTMPDIR':FTMPDIR,                  \
                         'platform':self.platform, 'output_root':self.output_root, 'diagTable':diagTable, 'dataTable':dataTable, \
                         'fieldTable':fieldTable, 'gridSpec':gridSpec,  'initCond':initCond, 'diagTable2':diagTable2, 'dataTable2':dataTable2, \
                         'fieldTable2':fieldTable2, 'gridSpec2':gridSpec2,  'initCond2':initCond2,  'inputData':inputData,              \
                         'npes':self.nnodes*self.ncores*self.nthreads, 'monthslist':self.monthslist, 'dayslist':self.dayslist, 'dt_atmos':self.dt_atmos,  \
                         'compiler_info':self.compiler_info, 'other_time_stamp':other_time_stamp, 'calendar':calendar}
        fmspy.io.write_script(self.run_filename, os.path.join(".", "bin/run.template"), output_format, option='overwrite')

    #===================================================================================
    # execute the run script
    #===================================================================================
    def exec_run(self):
        # exec the run script
        print("\nEXECUTING \"%s\"" % self.run_filename)
        os.chdir(os.path.dirname(self.run_filename))
        os.system(self.queue_cmd + ' ./' + os.path.basename(self.run_filename))
        os.chdir(self.cwd)

    #===================================================================================
    # setup and execute the compile script
    #===================================================================================
    def compile(self):
        self.setup_mkmf()
        self.setup_compile()
        status = self.exec_compile()

    #===================================================================================
    # setup and execute the run script
    #===================================================================================
    def run(self):
        self.setup_run()
        self.exec_run()

    #===================================================================================
    # setup and execute the compile and run scripts
    #===================================================================================
    def setup(self):
        self.setup_mkmf()
        self.setup_compile()
        status = self.exec_compile()
        if(status == 0): 
            self.setup_run()
            self.exec_run()

    #===================================================================================
    # setup and execute the compile and run scripts, but the run script is not executed
    #===================================================================================
    def setup2(self):
        self.setup_mkmf()
        self.setup_compile()
        status = self.exec_compile()
        if(status == 0): 
            self.setup_run()
