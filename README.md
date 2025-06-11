File thesis.pdf contains my bachelor thesis.

Next, the underlying codes for producing tables with results have been uploaded. Dependencies of individual functions are depicted in the figure function_dependencies.png.

The main script for simulating is simulacni_studie.R. This is the only place where parameters are to be modified.

It is possible for the output from simulacni_studie.R to be processed by the script simulacni_studie_out.R. However, this method was designed primarily for simulating via computational cluster, where large amounts of data could have been produced. For single-computer usage there is no need to dig into this script.

Note: time measuring was tried to be implemented. However, during runs of simulations, this information showed to be uninteresting. Hence, the implementation in some parts of the code is missing.
