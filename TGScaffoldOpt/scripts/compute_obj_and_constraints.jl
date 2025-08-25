using TGScaffoldOpt
using ForwardDiff
using BenchmarkTools

n = 10
# Params
D = 1000;
kf = 20.0;
λ = 0.00;
ρ₀ = 0.05;
growth_dir = "inward";
Tmax = 22;
# Initial Radius
# myR = 100ones(n);
θ0 = collect(range(0, 2π, length=n+1))
pop!(θ0)
myR = initial_structure("square",100,θ0,n);

function r_to_output(myR)
  D = 1000
  akf = 0.5
  bkf = 0.1
  A = 0.00
  ρ₀ = 0.05
  growth_dir = "inward"
  Tmax = 22

  θ, R, ρ = TGScaffoldOpt.TG_PDE_Solver(D, akf, bkf, A, ρ₀, Tmax, growth_dir, myR)
  j = mean(R[end])
end

function r_to_p_to_Vol(myR)
  p = TGScaffoldOpt.polar_to_cartesian(myR,θ0)
  Vol(p)
end

function r_to_p_to_Per(myR)
  p = TGScaffoldOpt.polar_to_cartesian(myR,θ0)
  Per(p)
end

# @benchmark r_to_output($myR)
# @benchmark ForwardDiff.gradient($r_to_output, $myR)

j = r_to_output(myR)
c1 = r_to_p_to_Vol(myR)
c2 = r_to_p_to_Per(myR)

dj = ForwardDiff.gradient(r_to_output, myR)
dC1 = ForwardDiff.gradient(r_to_p_to_Vol, myR)
dC2 = ForwardDiff.gradient(r_to_p_to_Per, myR)
# dVol = dC[1, :]
# dPer = dC[2, :]

using NLopt


function J(x::Vector, grad::Vector)
  if length(grad) > 0
    _dj = ForwardDiff.gradient(r_to_output,x)
    copyto!(grad,_dj)
  end
  return r_to_output(x)
end
function _Vol(x::Vector, grad::Vector)
  if length(grad) > 0
    _dj = ForwardDiff.gradient(r_to_p_to_Vol,x)
    copyto!(grad,_dj)
  end
  return r_to_p_to_Vol(x)
end
function _Per(x::Vector, grad::Vector)
  if length(grad) > 0
    _dj = ForwardDiff.gradient(r_to_p_to_Per,x)
    copyto!(grad,_dj)
  end
  return r_to_p_to_Per(x)
end

opt = NLopt.Opt(:LD_MMA, n)
NLopt.nlopt_set_param(opt, "verbosity", 1)
NLopt.lower_bounds!(opt, 0)
NLopt.xtol_rel!(opt, 1e-4)
NLopt.min_objective!(opt, J)
NLopt.inequality_constraint!(opt, (x, ∇) -> (-_Vol(x, ∇) + r_to_p_to_Vol(myR))/10, 1e-8)
NLopt.inequality_constraint!(opt, (x, ∇) -> -_Per(x, ∇) + 2π*sqrt(r_to_p_to_Vol(myR)/π), 1e-8)
min_f, min_x, ret = NLopt.optimize(opt, myR)
