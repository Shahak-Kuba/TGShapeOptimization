module SolverTest

using TGScaffoldOpt
using ForwardDiff

myR = 100ones(100)
T = eltype(myR)
D = 1000
kf = 20.0
λ = 0.00
ρ₀ = 0.05
growth_dir = "inward"
Tmax = 22
akf = 3.2741*1e-6
bkf = 8.5728*1e-5

@time θ, R, ρ = TGScaffoldOpt.TG_PDE_Solver(D, akf, bkf, λ, ρ₀, Tmax, growth_dir, myR)



function r_to_output(myR)
  T = eltype(myR)
  D = 1000
  kf = 20.0
  A = 0.00
  ρ₀ = 0.05
  growth_dir = "inward"
  Tmax = 22

  θ, R, ρ = TGScaffoldOpt.TG_PDE_Solver(T, D, kf, A, ρ₀, Tmax, growth_dir, myR)
  ρ[end]
end

@time ForwardDiff.gradient(r_to_output, myR)

@test true # It ran...

end