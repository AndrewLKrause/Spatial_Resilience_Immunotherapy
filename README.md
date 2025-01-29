TODO README

To run the files for continuation, you first have to download pde2path (REFERENCE HERE) and load it into Matlab.

In the folder 'pde2path continuation', there are two folders. Both of them work in the same way so the process will be explained for the folder called 'short domain' only.

Thereare two iportant steady states that are important for the parameter values chosen (the ones in Figure 3): the cancer-free equilibrium and another one for which v > 0.

To run the code for the cancer-free steady state, you must set "dir = 'p0';" in line 4 of the file 'Cancerinit.m'. With this, you can run the file cmds.m in order to get the continuation for the cancer-free steady state (In reality, only the first part of the computation matters in this case, as no branch points were continued to produce Figure 3).

After doing that, the computation will present an error because the folder 'p' has not been created and the files to plot the bifurcation curves associated with the other steady state will not be available.

Therefore, you must go back to the file 'Cancerinit.m', set "dir = 'p';", and run the file 'cmds.m' again.
