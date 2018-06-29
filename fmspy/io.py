#===================================================================================
# io.py: module to handle i/o in the fmspy package
#===================================================================================
import os
import sys
import shutil
sys.path.append("./fmspy")
import f90nml

SRC_TYPES       = ['.F90', '.f90', '.F', '.h', '.H', '.inc', '.c']
DOC_TYPES       = ['.HTML', '.html', '.dat', '.inp', '.nc',  '.txt', '.md', '.ps', '.pdf', '.doc', '.docx','.py', 'README', 'readme']
DAT_TYPES       = ['tar.gz', '.tar', '.dat', '.inp', '.nc', 'inputData.csh']

#===================================================================================
# copy one file from full_src_name to full_dst_name
# return the number of files copied
#===================================================================================
def copyfile(full_src_name, full_dst_name, option=None):

    if (os.path.exists(full_dst_name)):
        if option is None:
            print("\"%s\" already exists!" % full_dst_name)
            return 0
        elif option == 'silent':
            return 0
        elif option == 'overwrite':
            print("overwriting \"%s\"!" % full_dst_name)
    else:
        if not os.path.exists(os.path.dirname(full_dst_name)):
            os.makedirs(os.path.dirname(full_dst_name))
        print("copying: \"%s\"" % full_src_name)

    shutil.copy2(full_src_name, full_dst_name)

    return 1

#===================================================================================
# copy all the files or files of given file_types in the src folder to the dst folder
# return the number of files copied
#===================================================================================
def copyfolder(src, dst, files_to_exclude=None, file_types_to_copy=None, option=None):

    counts = 0
    src_files = os.listdir(src)
    if files_to_exclude is not None:
        for file_name in files_to_exclude:
            src_files.remove(file_name)

    for file_name in src_files:
        if file_types_to_copy is None:
            full_src_name = os.path.join(src, file_name)
            full_dst_name = os.path.join(dst, file_name)
            counts += copyfile(full_src_name, full_dst_name, option=option)
        else:
            for file_type in file_types_to_copy:
                if file_name.endswith(file_type):
                    full_src_name = os.path.join(src, file_name)
                    full_dst_name = os.path.join(dst, file_name)
                    counts += copyfile(full_src_name, full_dst_name, option=option)

    return counts

#===================================================================================
# copy src_files in the src folder and associated documentation files
# of given file_types (e.g., .html) in the same folders to the dst folder
# return the number of files copied
#===================================================================================
def copyfolder2(src, src_files, dst, file_types_to_copy=None, option=None):

    counts = 0
    for file_name in src_files:
        full_src_name = os.path.join(src, file_name)
        full_dst_name = os.path.join(dst, file_name)
        if file_types_to_copy is not None:
            # copy associated files of given file_types
            if not os.path.exists(os.path.dirname(full_dst_name)):
                os.makedirs(os.path.dirname(full_dst_name))
                copyfolder(os.path.dirname(full_src_name), os.path.dirname(full_dst_name), \
                    file_types_to_copy=file_types_to_copy, option=option)

        counts += copyfile(full_src_name, full_dst_name, option=option)

    return counts

#===================================================================================
# copy the source code from src to dst using pathnames in src and associated 
# files of given file_type (i.e., DOCT_TYPES) in the same folders
# return the number of files copied
#===================================================================================
def copysrc(src, dst, option=None):

    # make src directory
    print("\nsetting up \"%s\"" % dst)
    if not os.path.exists(dst):
        os.makedirs(dst)

    # copy src files
    src_files = get_src_files(src, file_types=SRC_TYPES, root_to_skip="./my", root_to_include="", file_to_skip=[])
    counts = copyfolder2(src, src_files, dst, file_types_to_copy=DOC_TYPES, option=option)
    print("Total number of src files copied: \"%d\"\n" % counts)

    return counts

#===================================================================================
# create mkmf, compile, and run scripts with template and output_format
# write the output in the file script_name
#===================================================================================
def write_script(script_name, template, output_format, option=None):

    if os.path.exists(script_name):
        if option is None:
            print("\"%s\" already exists!" % script_name)
        elif option == "overwrite":
            print("overwriting \"%s\"!" % script_name)
    else:
        print("creating \"%s\"" % script_name)

    # read the generic template
    f = open(template, 'r')
    output = f.read()
    f.close()
    
    # write the script file with output_format
    #print(output % output_format)
    f = open(script_name, 'w')
    f.write(output % output_format)
    f.close()

    # change file permission
    os.chmod(script_name, 0o775)


#===================================================================================
# generate a list of pathnames from the src directory
# return a list of src_files
#===================================================================================
def get_src_files(src, file_types=SRC_TYPES, root_to_skip="./my", root_to_include="", file_to_skip=[]):

    cwd = os.getcwd()
    os.chdir(src)

    src_files = []
    for root, dirs, files in os.walk("./"):
        if not root.startswith(root_to_skip):
            if root.startswith(root_to_include):
                for file in files:
                    for file_type in file_types:
                        if file.endswith(file_type):
                            full_file = os.path.join(root,file).replace("./", "")
                            if not full_file in file_to_skip: 
                                src_files.append(full_file)

    src_files.sort()
    os.chdir(cwd)

    return src_files

#===================================================================================
# src files in src1 but not in src2
#===================================================================================
def src_files_diff(src1, src2, file_types=SRC_TYPES, root_to_skip1="./my", root_to_skip2="./my", \
    root_to_include1="", root_to_include2="", file_to_skip1=[], file_to_skip2=[]):

    src_files1  = get_src_files(src1, file_types=file_types, root_to_skip=root_to_skip1, root_to_include=root_to_include1, file_to_skip=file_to_skip1)
    src_files2  = get_src_files(src2, file_types=file_types, root_to_skip=root_to_skip2, root_to_include=root_to_include2, file_to_skip=file_to_skip2)
    src_diff    = [file for file in src_files1 if file not in src_files2]

    return src_diff

#===================================================================================
# write the path_names file using src_files
#===================================================================================
def write_src_files(dst, src_files):

    cwd      = os.getcwd()
    full_dst = os.path.abspath(dst)
    os.chdir(os.path.dirname(full_dst))

    f = open(os.path.basename(full_dst), 'w')
    for filename in src_files:
        f.write("%s\n" % filename)
    f.close()

    os.chdir(cwd)
    print("Total number of src files written in path_names: \"%d\"\n" % len(src_files))

#===================================================================================
# copy diag_table, data_table, field_table, and namelists in one experiment from src to dst
#===================================================================================
def copy_exp_setup(src, dst, namelist_patch={}, option='overwrite'):

    if (not os.path.exists(dst)):
        os.makedirs(dst)

    if os.path.exists(src+'/diag_table'):
        copyfile(src+'/diag_table', dst+'/diag_table', option=option)
        
    if os.path.exists(src+'/data_table'):
        copyfile(src+'/data_table', dst+'/data_table', option=option)
        
    if os.path.exists(src+'/field_table'):
        copyfile(src+'/field_table', dst+'/field_table', option=option)
        
    # copy the namelist file using the f90nml package 
    full_src_name = os.path.join(src, 'namelists')
    full_dst_name = os.path.join(dst, 'namelists')
    if namelist_patch:
        print("Namelists in %s are overriden: %s" % (full_dst_name, namelist_patch))
        if os.path.exists(full_dst_name):
            os.remove(full_dst_name)
        f90nml.patch(full_src_name, namelist_patch, full_dst_name)
    else:
        copyfile(full_src_name, full_dst_name, option=option)

#===================================================================================
# length of integration using monthslist, dayslist, secondslist
# default: monthslist, dayslist, secondslist = "(0)", "(1)", "(0)"
#===================================================================================
def months_days_seconds_list(monthslist, dayslist, secondslist):
    def str2tuple(s):
        if (s == None):
            return([])
        else:
            si = s.find('(')
            se = s.find(')')
            ss = s[si+1:se].strip()
        return(ss.split())

    mlen = len(str2tuple(monthslist))
    dlen = len(str2tuple(dayslist))
    slen = len(str2tuple(secondslist))
    ilen = max(mlen, dlen, slen)

    if (ilen == 0):
        #default: integrating the model by 1 day
        monthslist  = "(0)"
        dayslist    = "(1)"
        secondslist = "(0)"
    else:
        if(mlen < ilen):
            monthslist  = "( " + "0 "*ilen + ")"
        if(dlen < ilen):
            dayslist    = "( " + "0 "*ilen + ")"
        if(slen < ilen):
            secondslist = "( " + "0 "*ilen + ")"

    return monthslist, dayslist, secondslist
