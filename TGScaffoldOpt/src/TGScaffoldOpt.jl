module TGScaffoldOpt
using LinearAlgebra
using ForwardDiff
using CairoMakie

include("PDE_Solver.jl")
export TG_PDE_Solver

include("ConstraintsCalcs.jl")
export Vol, Per

include("Utils.jl")
export initial_structure, plot_iteration_TG_sol

end
