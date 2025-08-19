module TGScaffoldOpt
using LinearAlgebra
using ForwardDiff
using CairoMakie

include("PDE_Solver.jl")
include("ConstraintsCalcs.jl")
include("Visualisation.jl")

export TG_PDE_Solver

end # module TGScaffoldOpt
