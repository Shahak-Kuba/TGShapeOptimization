using TGScaffoldOpt
using ForwardDiff
using BenchmarkTools
using TGScaffoldOpt: polar_to_cartesian, Vol, Per

myR = 100ones(100);
D = 1000;
kf = 20.0;
λ = 0.00;
ρ₀ = 0.05;
growth_dir = "inward";
Tmax = 22;

function r_to_output(myR)
  D = 1000
  kf = 20.0
  A = 0.00
  ρ₀ = 0.05
  growth_dir = "inward"
  Tmax = 22

  θ, R, ρ = TGScaffoldOpt.TG_PDE_Solver(D, kf, A, ρ₀, Tmax, growth_dir, myR)
  j = ρ[end]
end

function r_to_constraints(myR)
  m = size(myR, 1) + 1
  θ = collect(range(0, 2π, length=m))
  pop!(θ)
  p = polar_to_cartesian(myR,θ)
  [Vol(p), Per(p)]
end

# @benchmark r_to_output($myR)
# @benchmark ForwardDiff.gradient($r_to_output, $myR)

j = r_to_output(myR)
c = r_to_constraints(myR)

dj = ForwardDiff.gradient(r_to_output, myR)
dC = ForwardDiff.jacobian(r_to_constraints, myR)
dVol = dC[1, :]
dPer = dC[2, :]