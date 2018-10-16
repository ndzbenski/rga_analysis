# rga_analysis

For use: clone the repository. Make sure you have Groovy installed, then do
> <clara-directory>/installation/plugins/clas12/bin/run-groovy <groovy-file> <hipo-file(s)> [beam energy]

where <clara-directory> is where the working version of clara is installed, <groovy-file> can be any of the groovy analysis files from this directory, <hipo-file(s)> is an individual hipo file or a file containing a list of hipo files for analysis, and finally the beam energy of the hipo files.

For example, if you want to run an uncut analysis on a set of cooked files, do
> <clara-directory>/installation/plugins/clas12/bin/run-groovy uncut_analysis.groovy cooked.files 10.6
