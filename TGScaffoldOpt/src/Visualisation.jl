using CairoMakie


function PlotIterationTGSol(R,θ,ρ,sol_indecies)
    fig = Figure(size=(600, 600))
    ax = Axis(fig[1, 1], xlabel="X", ylabel="Y")

    for idx in sol_indecies
        Rₙ = [R[idx, :]; R[idx, 1]]  # Close the loop by repeating the first point
        ρₙ = [ρ[idx, :]; ρ[idx, 1]]  # Close the loop by repeating the first point
        #Rₙ = R[idx, :]  
        #ρₙ = ρ[idx, :]
        coords = polar_to_cartesian(Rₙ, [θ; 2π])
        lines!(ax, coords[1, :], coords[2, :], color=ρₙ, colormap=:viridis, linewidth=3, colorrange=[0,0.08])
    end

    return fig
end

