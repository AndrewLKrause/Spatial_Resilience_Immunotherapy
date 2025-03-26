This GitHub collects code for linear stability analyses, direct 1D and 2D simulations, and 1D numerical continuation for the paper, "Pattern Formation as a Resilience Mechanism in Cancer Immunotherapy." We also refer to the VisualPDE page on the [immunotherapy model](https://visualpde.com/mathematical-biology/immunotherapy-model) for an interactive implementation of the PDE models.

## Linear Stability Analyses
The codes `LinStabODE.m` and `TuringStability.m` respectively perform the linear stability analysis described in the paper for the ODE and PDE models respectively, including producing the figures. The file `Symbolic.m' was used to check some of the analytical calculations in these analyses.

## Direct simulations
The code implementing the simulations is `Solve_Model.m` where you set the dimension, `dim`, and number of gridpoints, `m`, to use in the method-of-lines discretization. Other model parameters are kept in the file `Init_Parameters.m`.

## Numerical Continuation via pde2path
To run the files for continuation, you first have to download pde2path (see [https://pde2path.uol.de/index.html](https://pde2path.uol.de/index.html)) and load it into the current path in Matlab.

In the folder `pde2path continuation`, there are two folders. Both of them work in the same way so the process will be explained for the folder called `short domain` only.

We focus on continuing along two kinds of spatially homogeneous equilibria: the cancer-free equilibrium and the cancer coexistence equilibrium for which $v > 0$.

To run the code for the cancer-free steady state, you must set `dir = 'p0';` in line 4 of the file `Cancerinit.m`. With this, you can run the file `cmds.m` in order to obtain the continuation for the cancer-free steady state (In reality, only the first part of the computation matters in this case, as no branch points were continued to produce the Figures).

You must then go back to the file `Cancerinit.m`, set `dir = 'p';`, and run the file `cmds.m` again.
