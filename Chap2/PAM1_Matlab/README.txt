Instructions on running the solver. This was run with Matlab 2015b, and depends
on very few built-in functions --> Obsolescence unlikely.

1) For quick simulation, type <RUNFILE> into Matlab console and wait for results

Explanation of files within folder.

stoichBalance.ods: balancing of the stoichiometry to ensure mass balances were
met over C, COD, N and P. The results of this balance were input into the
PAM1.m function file.

dynamicInfluent.csv: influent characteristics to simulate a PANMBR continuously
over a period of ~ 600 days.

pHsolve.m: Function to solve the pH of the system. Is a dependency of PAM1.m

balances.m: Run after the simulation is completed to check that the mass balance
checks out. This is just a sanity check. 

state... : results file for both influent and effluent.

PAM1.m: The function file. Contains all important information of the PPB system.
This is the file to be modified if improvements are to be made.

RUNFILE.m: Glues everything together. To change simulation conditions such as
run time, integrator, graphing preferences, etc, modify this file.
