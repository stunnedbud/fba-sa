# This file requires the package COBRA. It can be installed with the command Pkg.add("COBRA")
# COBRA requires a Julia version of 0.5 or above.
using COBRA

# include the solver configuration file
include("$(Pkg.dir("COBRA"))/config/solverCfg.jl")

# COBRA requires a solver. In this file the solver Clp is used (Pkg.add("Clp")).
# change the COBRA solver to Clp
solver = changeCobraSolver("Clp", [])

# Load the stoichiometric matrix S from a MATLAB structure named model in the specified .mat file
#model = loadModel("sources/ecoli_core_model.mat", "S", "model");

# TO DO: get the script running using this file, the regulated model ie. with all constraints?
model = loadModel("sources/modelReg.mat", "S", "modelReg");

# set the reaction list (only one reaction)
rxnsList = 13

# select the reaction optimization mode
#  0: only minimization
#  1: only maximization
#  2: maximization and minimization
rxnsOptMode = 1

# launch the distributedFBA process with only 1 reaction on 1 worker
minFlux, maxFlux, optSol, fbaSol, fvamin, fvamax = distributedFBA(model, solver, nWorkers=1, rxnsList=rxnsList, rxnsOptMode=rxnsOptMode);

