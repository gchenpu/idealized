#===================================================================================
# main program
#===================================================================================
import os
import fmspy.io

#===================================================================================
model_name     = "f90"
version        = "Tpdf"
runscript_name = "Tpdf"

#model_name      = "hs_with_clouds"
#version         = "default"
#runscript_name  = "hs_with_clouds"

#===================================================================================
if(not os.path.exists("../data/release/"+model_name+"/exp/"+model_name+"/"+version)):
    os.makedirs("../data/release/"+model_name+"/exp/"+model_name+"/"+version)

# copy files
os.system("cp -r docs fmspy netcdf postprocessing ../data/release/"+model_name)
os.system("cp runpy/run_"+runscript_name+".py ../data/release/"+model_name)
os.system("cp -r ./bin "+" ../data/release/"+model_name)
os.system("cp -r ./exp/"+model_name+"/input_data"+" ../data/release/"+model_name+"/exp/"+ model_name+"/")
os.system("cp -r ./exp/"+model_name+"/"+ version +" ../data/release/"+model_name+"/exp/"+ model_name+"/")

# copy src code
src        = "./src"
dst        = "../data/release/"+model_name+"/src"
pathnames  = "./exp/"+model_name+"/"+version+"/src/path_names"

f = open(pathnames, 'r')
src_files = f.read().split()
f.close()

counts = fmspy.io.copyfolder2(src, src_files, dst, file_types_to_copy=fmspy.io.DOC_TYPES, option='silent')
print("Total number of src files copied: \"%d\"\n" % counts)

# create tar file
os.chdir("../data/release/")
os.system("tar -cvf "+model_name+".tar "+model_name)
