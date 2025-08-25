## Initial structure
function initial_structure(type, r₀, θ, M)
  T = eltype(θ)
  r = zeros(T, M)

  if type == "circle"
    r .= r₀
  elseif type == "square"
    R = T(50.0)
    for i in 1:M
      if θ[i] < π / 2
        r[i] = min(R / cos(θ[i]), R / sin(θ[i]))
      else
        r[i] = min(R / cos(mod(θ[i], π / 2)), R / sin(mod(θ[i], π / 2)))
      end
    end
  elseif type == "hex"
    R = sqrt((2 / (3 * sqrt(3))) * π * (r₀^2))
    for i in 1:M
      if θ[i] >= 0 && θ[i] < π / 3
        r[i] = R / (cos(θ[i]) + 1 / sqrt(3) * sin(θ[i]))
      elseif θ[i] >= π / 3 && θ[i] < π / 2
        r[i] = R * sqrt(3) / (2 * sin(θ[i]))
      elseif θ[i] >= π / 2 && θ[i] < 2 * π / 3
        r[i] = R * sqrt(3) / (2 * sin(θ[i]))
      elseif θ[i] >= 2 * π / 3 && θ[i] < π
        r[i] = R / (-cos(θ[i]) + 1 / sqrt(3) * sin(θ[i]))
      elseif θ[i] >= π && θ[i] < 4 * π / 3
        r[i] = -R / (cos(θ[i]) + 1 / sqrt(3) * sin(θ[i]))
      elseif θ[i] >= 4 * π / 3 && θ[i] < 3 * π / 2
        r[i] = -R * sqrt(3) / (2 * sin(θ[i]))
      elseif θ[i] >= 3 * π / 2 && θ[i] < 5 * π / 3
        r[i] = -R * sqrt(3) / (2 * sin(θ[i]))
      else
        r[i] = R / (cos(θ[i]) - 1 / sqrt(3) * sin(θ[i]))
      end
    end
  end
  return r
end

## Visualisation
function plot_iteration_TG_sol(R, θ, ρ, sol_indecies;
      cmap=:viridis, lw=3, ClrRange=[0, 0.08])
  fig = Figure(size = (600,600))
  ax = Axis(fig[1, 1], xlabel="X", ylabel="Y")

  for idx in sol_indecies
    Rₙ = [R[idx, :]; R[idx, 1]]  # Close the loop by repeating the first point
    ρₙ = [ρ[idx, :]; ρ[idx, 1]]  # Close the loop by repeating the first point
    #Rₙ = R[idx, :]
    #ρₙ = ρ[idx, :]
    coords = polar_to_cartesian(Rₙ, [θ; 2π])
    lines!(ax, coords[1, :], coords[2, :], ρₙ, colormap = cmap, linewidth = lw, colorrange = ClrRange)
  end

  return fig
end