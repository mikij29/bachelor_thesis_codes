# Robust Linear Regression Simulation Interface 2.0
This is a set of R-scripts for simulating linear regression data and applying least weighted squares (LWS), among many other (not only robust) estimators on them. Among implemented features there are various metrics of errors regarding the accuracy aspect of estimators. Next, two methods for computing their variance are used - sandwich method and nonparametric bootstrap. 

## Files
File 'thesis.pdf' contains my bachelor thesis, for which the codes have been created primarily. Old version of these codes was handed in as an attachment of that thesis, but contained several mistakes and inconsistencies. Newer, but still somehow problematic versions were saved here as a backup.

### R-Scripts
Dependencies of individual functions (which are usually located in their own files) are depicted in the figure 'function_dependencies.png'. Note: functions from 'load_names.R' (the only one not contained in the graph) are auxiliary and used all throughout the code esp. when calculating dimensions is needed.

The main script for simulating is 'simulacni_studie.R'. This is the only place where parameters are to be modified.

#### Input 
On the input, there are plenty of parameters to modify. Except of how to generate data (and how to contaminate them) and which estimators to use, the user might set for example their own number of iterations used for initial estimation in data-dependent weighted least squares or data-dependent LWS procedures.

#### Output
Output of all the procedures might be chosen either as a console printed table, or as a latex file containing the table as well as some of its formatting. Abbreviations 'sand', 'boot' and 'jest' stand for sandwich method, bootstrap and 'just estimate' (estimation without corresponding estimate of variance) respectively.

Note: time measuring was tried to be implemented. However, during runs of simulations, this information showed as uninteresting. Hence, the implementation in some parts of the code is still missing.

### Files for simulating via cluster
Most of the simulations were performed via running a shell script in form as in qsim1.sh. 'Ranks' file is as simple as a list of numbers, as shown in ranks200.txt. Note: both files need unix-friendly End Of Line characters to work properly. Hence, after creating files like the ones above somehow in Windows, I needed to convert them by means of command ":set fileformat=unix" in Vim editor.

## Requirements
You need to install R and its library 'MASS'. The library 'xtable' is needed for producing latex tables in the output.

More requirements are needed when simulating via a computational cluster. I used the one with homepage https://cluster.karlin.mff.cuni.cz/ and would like to thank the administrators for helping me with my beginnings.

## Author
Michal Jelínek, Faculty of Mathematics and Physics - Department of Probability and Mathematical Statistics, Charles University, Prague.

Do not hesitate to contact me (jelinek@karlin.mff.cuni.cz) with any questions or issues. Eventually, I will be happy to adjust the codes based on your preferences. 
