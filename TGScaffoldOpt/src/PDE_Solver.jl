# Semi-implicit finite difference solver for governing PDE

# TG PDE is given by: ∂ρ/∂t = D(∂²ρ/∂x²) + ρvκ + λρ 
# where v = k_f ρ, D is diffusion, κ is curvature, λ is a population growth/decay term.

using LinearAlgebra
# In-place upwind selector (Dual-safe)
function aₘ!(a⁺::AbstractVector{T}, a⁻::AbstractVector{T},
             Rᵢ₊₁::AbstractVector{T}, Rᵢ₋₁::AbstractVector{T}) where {T}
    oneT, zeroT = oneunit(T), zero(T)
    @inbounds @simd for i in eachindex(a⁺)
        ap = ifelse(Rᵢ₊₁[i] - Rᵢ₋₁[i] > zeroT, oneT, zeroT)
        a⁺[i] = ap
        a⁻[i] = oneT - ap
    end
    return a⁺, a⁻
end

"""
    tridiag_solve!(x, dl, d, du, b, cwork, dwork)

Solve A*x = b where A has lower diag `dl` (n-1), main diag `d` (n),
upper diag `du` (n-1). Uses work arrays `cwork`, `dwork` (length n).
Writes solution into `x`. Dual-safe.
"""
function tridiag_solve!(x, dl, d, du, b, cwork, dwork)
    n = length(d); oneT = oneunit(eltype(d))
    @inbounds begin
        invd = oneT / d[1]
        cwork[1] = du[1] * invd
        dwork[1] = b[1]  * invd
        for i in 2:n-1
            denom  = d[i] - dl[i-1] * cwork[i-1]
            invden = oneT / denom
            cwork[i] = du[i] * invden
            dwork[i] = (b[i] - dl[i-1] * dwork[i-1]) * invden
        end
        denom = d[n] - dl[n-1] * cwork[n-1]
        dwork[n] = (b[n] - dl[n-1] * dwork[n-1]) / denom

        x[n] = dwork[n]
        for i in (n-1):-1:1
            x[i] = dwork[i] - cwork[i] * x[i+1]
        end
    end
    return x
end

"""
    cyclic_tridiag_solve!(x, dl, d, du, α, β, b, cwork, dwork, z1, zN, y)

Solve (A + corners) * x = b where A has bands `dl` (n-1), `d` (n), `du` (n-1),
and corner couplings A[1,n] = α, A[n,1] = β. Dual-safe. No allocations.
`cwork`, `dwork`, `z1`, `zN`, `y` are length-n work vectors.
"""
function cyclic_tridiag_solve!(x, dl, d, du, α, β, b, cwork, dwork, z1, zN, y)
    T = eltype(d); oneT, zeroT = oneunit(T), zero(T); n = length(d)

    tridiag_solve!(y,  dl, d, du, b,  cwork, dwork)   # y = A \ b
    fill!(z1, zeroT); z1[1]  = oneT
    tridiag_solve!(z1, dl, d, du, z1, cwork, dwork)   # z1 = A \ e1
    fill!(zN, zeroT); zN[n]  = oneT
    tridiag_solve!(zN, dl, d, du, zN, cwork, dwork)   # zN = A \ eN

    s11 = oneT + α*z1[n]; s12 =      β*zN[n]
    s21 =      α*z1[1];   s22 = oneT + β*zN[1]
    r1, r2 = y[n], y[1]
    detS = s11*s22 - s12*s21
    c1 = ( s22*r1 - s12*r2) / detS
    c2 = (-s21*r1 + s11*r2) / detS

    @inbounds @simd for i in 1:n
        x[i] = y[i] - α*z1[i]*c1 - β*zN[i]*c2
    end
    return x
end

# --------------------------
# Main solver
# --------------------------

function TG_PDE_Solver(T, D, kf, A, ρ₀, Tmax, growth_dir, myR)
    # grid (keep scalars in T)
    T₀  = zero(T)
    N   = 10000
    Δt  = (T(Tmax) - T₀) / T(N)
    m   = 101
    Δθ  = T(2π) / T(m)
    θ   = collect(range(T(0), T(2π), length=m)); pop!(θ)
    M   = m - 1
    invΔθ  = oneunit(T) / Δθ
    invΔθ2 = invΔθ * invΔθ

    # state arrays
    R   = zeros(T, N+1, M)
    Rθ  = zeros(T, N+1, M)
    Rθθ = zeros(T, N+1, M)
    ρ   = zeros(T, N+1, M)
    ρθ  = zeros(T, N+1, M)
    κ   = zeros(T, N+1, M)
    Φ   = zeros(T, N+1, M)
    λ   = zeros(T, N+1, M)
    ψ   = zeros(T, N+1, M)

    # initial conditions
    @views copyto!(R[1, :], myR)
    @views ρ[1, :] .= ρ₀

    S = growth_dir == "inward" ? oneunit(T) : -oneunit(T)

    # scratch (reused)
    a⁺   = Vector{T}(undef, M)
    a⁻   = Vector{T}(undef, M)
    Rᵢ   = Vector{T}(undef, M)
    ρᵢ   = Vector{T}(undef, M)
    Rᵢ₊₁ = Vector{T}(undef, M)
    Rᵢ₋₁ = Vector{T}(undef, M)
    ρᵢ₊₁ = Vector{T}(undef, M)
    ρᵢ₋₁ = Vector{T}(undef, M)

    # tridiagonal bands & work
    d₀   = Vector{T}(undef, M)
    d₁   = Vector{T}(undef, M-1)
    d₋₁  = Vector{T}(undef, M-1)
    cwrk = Vector{T}(undef, M)
    dwrk = Vector{T}(undef, M)
    z1   = Vector{T}(undef, M)
    zN   = Vector{T}(undef, M)
    ybuf = Vector{T}(undef, M)
    xsol = Vector{T}(undef, M)

    @views @inbounds for n in 1:N
        Rn    = R[n,   :];  Rnp1 = R[n+1, :]
        ρn    = ρ[n,   :];  ρnp1 = ρ[n+1, :]
        Rθn   = Rθ[n,  :];  Rθθn = Rθθ[n, :]
        ρθn   = ρθ[n,  :]
        κn    = κ[n,   :]
        Φn    = Φ[n,   :]
        λn    = λ[n,   :]
        ψn    = ψ[n,   :]

        copyto!(Rᵢ, Rn)
        copyto!(ρᵢ, ρn)

        # neighbors (no circshift allocs)
        @inbounds for i in 1:M
            ip = (i == M) ? 1 : i+1
            im = (i == 1) ? M : i-1
            Rᵢ₊₁[i] = Rᵢ[ip]; Rᵢ₋₁[i] = Rᵢ[im]
            ρᵢ₊₁[i] = ρᵢ[ip]; ρᵢ₋₁[i] = ρᵢ[im]
        end

        # upwind flags
        aₘ!(a⁺, a⁻, Rᵢ₊₁, Rᵢ₋₁)

        # derivatives + geometry + coefficients (fused)
        @inbounds for i in 1:M
            Rθi  = (Rᵢ[i]*(a⁺[i]-a⁻[i]) - Rᵢ₋₁[i]*a⁺[i] + Rᵢ₊₁[i]*a⁻[i]) * invΔθ
            Rθθi = (Rᵢ₊₁[i] - (T(2)*Rᵢ[i]) + Rᵢ₋₁[i]) * invΔθ2
            ρθi  = (ρᵢ[i]*(a⁺[i]-a⁻[i]) - ρᵢ₋₁[i]*a⁺[i] + ρᵢ₊₁[i]*a⁻[i]) * invΔθ

            Rθn[i]  = Rθi
            Rθθn[i] = Rθθi
            ρθn[i]  = ρθi

            R2   = Rᵢ[i]*Rᵢ[i]
            Rθ2  = Rθi*Rθi
            denom = (R2 + Rθ2)^(T(1.5))
            κn[i] = -(R2 + T(2)*Rθ2 - Rᵢ[i]*Rθθi) / denom

            hypotR = sqrt(R2 + Rθ2)
            Φn[i] = ρᵢ[i] - S*Δt*κn[i]*kf*(ρᵢ[i]*ρᵢ[i]) -
                    (S*Δt*kf*ρᵢ[i]*Rθi*ρθi) / (Rᵢ[i]*hypotR) - A*Δt*ρᵢ[i]

            denomλ = R2 + Rθ2
            λn[i]  = (D*Δt*invΔθ2) / denomλ
            ψn[i]  = (D*Δt*Rθi*(Rᵢ[i] + Rθθi)) / (Δθ * (denomλ*denomλ))
        end

        # explicit update for R
        @inbounds @simd for i in 1:M
            Rnp1[i] = Rᵢ[i] - (S*Δt*kf*ρᵢ[i]/Rᵢ[i]) * sqrt(Rᵢ[i]*Rᵢ[i] + Rθn[i]*Rθn[i])
        end

        # build tridiagonal bands
        @inbounds for i in 1:M
            d₀[i] = oneunit(T) + T(2)*λn[i] + ψn[i]*(a⁺[i] - a⁻[i])
            if i < M
                d₁[i]  = -λn[i]   + ψn[i]*a⁻[i]
                d₋₁[i] = -λn[i+1] - ψn[i+1]*a⁺[i+1]
            end
        end
        α = -λn[1] - ψn[1]*a⁺[1]     # A[1,M]
        β = -λn[M] + ψn[M]*a⁻[M]     # A[M,1]

        # ----- cyclic tridiagonal solve  -----
        cyclic_tridiag_solve!(xsol, d₋₁, d₀, d₁, α, β, Φn, cwrk, dwrk, z1, zN, ybuf)
        copyto!(ρnp1, xsol)
    end

    return θ, R, ρ
end