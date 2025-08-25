using TGScaffoldOpt
using CairoMakie
using ForwardDiff
using BenchmarkTools
using TGScaffoldOpt: polar_to_cartesian, Vol, Per, PlotIterationTGSol

n = 400
# Params
D = 30;
kf = 20.0;
λ = 0.00;
ρ₀ = 0.05;
growth_dir = "inward";
Tmax = 23;
# Initial Radius
#myR = 100ones(n);
θ0 = collect(range(0, 2π, length=n+1))
pop!(θ0)
myR = r_for_InitialBoundary("square",100,θ0,n);

@time θ, R, ρ = TGScaffoldOpt.TG_PDE_Solver(D, kf, λ, ρ₀, Tmax, growth_dir, myR)

akf = 0.5#3.2741*1e-6
bkf = 0.1#8.5728*1e-5

@time θ, R, ρ = TGScaffoldOpt.TG_PDE_Solver(D, akf, bkf, λ, ρ₀, Tmax, growth_dir, myR)

f = TGScaffoldOpt.plot_iteration_TG_sol(R,θ,ρ,[1, 2001, 4001, 6001, 8001, 10001])


