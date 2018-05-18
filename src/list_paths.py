#===================================================================================
# create the path_names file for the src directory
#===================================================================================
import sys
import os
sys.path.append("..")
sys.path.append("../fmspy")
import fmspy.io

file_types     = fmspy.io.SRC_TYPES 
root_to_skip   = "./my"
file_to_skip   = [] 

src_files = fmspy.io.get_src_files("./", file_types=file_types, root_to_skip=root_to_skip, file_to_skip=file_to_skip)
fmspy.io.write_src_files("path_names", src_files)
