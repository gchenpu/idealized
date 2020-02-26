=========  How to run the model and modify the model code =========

=== basic info
idealized/exp/hs/default:  default experiment configurations and src code
idealized/src:             src code shared by different models; idealize/src/my will be overwritten by the src code in idealized/exp/hs/default/src
idealized/run_fms.py:      python script to compile and run the model
idealized/workdir:         the folder where temporatory files are saved


=== run the model under a new version name "test"
1. find the right model (e.g., hs)
   cd into exp/hs and create a new model version name "test"
   cp run_ctl and src in exp/hs/default/ to your new folder "test"
   run the python script in the root directory and specify your new version name [i.e., version='test']
   add the new version name to .gitignore unless you are ready to update it to the repository
   
=== modify the code   
2. add or modify your own scripts in the folder "test/src"; 
   this is the permenant directory to save changes to the model code, which will overwite idealized/src/my
   add a new experiment by creating test/run_test and copy everything in test/run_ctl to the folder
   *document your updates to the model by README.md
