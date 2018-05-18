# To retrieve a local copy of the remote server
git clone https://github.com/gchenpu/idealized.git

# suggested git gui tools
type "git gui"; this is a tool for managing commits
type "gitk"; this is a graphical history viewer
# see more at https://nathanj.github.io/gitguide/tour.html

# command line
# To update the code to the remote repository in a sequence
git add *
git commit -m ""
git push origin master

# Individual steps
# To propose changes and add them to the INDEX
git add filename
git add *

# To actually commit these changes to the HEAD
git commit -m "Commit message"

# To send those changes to your remote repository, execute (username & password will be prompted)
git push origin master

# To connect your repository to a remote server and add it
git remote add origin https://github.com/gchenpu/idealized.git

# See more at http://rogerdudler.github.io/git-guide/

