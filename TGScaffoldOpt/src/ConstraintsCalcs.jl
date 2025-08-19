"""
    polar_to_cartesian(R, θ)

Convert a row vector of radii `R` and a vector of angles `θ` into a matrix of 2D Cartesian coordinates.
Each column of the output matrix represents a position vector [x; y].

# Arguments
- `R`: A 1D array or row vector of radii.
- `θ`: A 1D array or vector of angles (in radians).

# Returns
A matrix where each column is a 2D position vector [x; y].
"""
function polar_to_cartesian(R::Vector{T}, θ::Vector{T}) where T
    n = length(R)
    coords = zeros(T,2,n)
    for i in 1:n
        coords[1, i] = R[i] * cos(θ[i])
        coords[2, i] = R[i] * sin(θ[i])
    end
    return coords
end

"""
    Ω(p)

Calculate the area of a polygon defined by points in `p`.

This function computes the area using the shoelace formula. The polygon is defined by a set of points `p`, and the function iterates through these points to calculate the area.

# Arguments
- `p`: A matrix where each column represents a point of the polygon in 2D space.

# Returns
The absolute area of the polygon.
"""
function Ω(p)
    A = 0
    for ii in axes(p,2)
        if ii == size(p,2)
            A += (p[1,ii]*p[2,1] -  p[2,ii]*p[1,1])
        else
            A += (p[1,ii]*p[2,ii+1] -  p[2,ii]*p[1,ii+1])
        end
    end
    return abs(A)/2;
end


"""
    P(p)

Calculate the perimeter of a polygon defined by points in `p`.

# Arguments
- `p`: A matrix where each column represents a point of the polygon in 2D space.

# Returns
The perimeter (sum of edge lengths) of the polygon.
"""
function P(p)
    n = size(p, 2)
    perim = 0.0
    for i in 1:n
        j = i == n ? 1 : i + 1
        perim += norm(p[:, i] - p[:, j])
    end
    return perim
end
