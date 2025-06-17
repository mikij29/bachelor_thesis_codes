File thesis.pdf contains my bachelor thesis.

R-SCRIPTS
Next, the underlying codes for producing tables with results have been uploaded. Dependencies of individual functions are depicted in the figure function_dependencies.png. Note: functions from load_names.R (the only one not contained in the graph) are auxiliary and used throughout the code esp. when calculating dimensions is needed.

The main script for simulating is simulacni_studie.R. This is the only place where parameters are to be modified.

Note: time measuring was tried to be implemented. However, during runs of simulations, this information showed as uninteresting. Hence, the implementation in some parts of the code is missing.

FILES FOR SIMULATING VIA CLUSTER
Most of the simulations were performed via running a shell script in form as in qsim1.sh. 'Ranks' file is as simple as a list of numbers, as shown in ranks200.txt. Note: both files need unix-friendly End Of Line characters to work properly. Hence, after creating files like the ones above somehow in Windows, I needed to convert them by means of command ":set fileformat=unix" in Vim editor.
