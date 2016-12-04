Code in R Creates Graphs, samples data and runs the FCI algorithm. Also, creates marginal data sets and runs FCI on them.

Files "fci.R" and "udag2pag.R" in folder codeFCI are copies from the pcalg R package, but we added the "aggressively prevent cycles" option. Thus, if one wants to have the option, they need to source these files after loading the package.

!!CAREFUL!! In the case of marginals, the outputs refer to the randomly ordered variables. Thus, if the vector "vars" is: vars=[3,7,...], it means that the first column/row of the adjacency matrix, separating sets etc. correspond to the 3rd variable of the original data (the one containing all variables). Re-ordering is done after the data analysis (matlab code).
