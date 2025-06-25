module solOZ

export mix_solve

include("legendreExpansion.jl")
using .legendreExpansion

function simpleIter(x::Vector{Float64}, N::Int64, gamma::Vector{Float64}, u::Vector{Float64}, bridge::Vector{Float64}, legendrePol::Vector{Vector{Float64}})
    # Using closure 
    c_dir = exp.(-1.0 .* u .+ gamma .+ bridge) .- gamma .- 1.0
    
    # Calculating Legendre coefficients for c_dir
    c_l = leg_coeff_arr(x, c_dir, legendrePol)
 
    l = collect(0:length(c_l)-1)
    
    # Applying OZ equation to obtain gamma_l
    denominator = 2.0 .* l .+ 1.0
    fraction = (N ./ denominator) .* c_l
    adjusted_c_l = c_l ./ (1.0 .- fraction)
    
    # Final expression
    gamma_l = adjusted_c_l .- c_l
    gamma_new = series_expansion(gamma_l, legendrePol)
    
    return gamma_new
end

"""
This function uses the simpleIter function to iteratively solve the OZ equation using a mixing parameter α.

"""
function mix_solve(α, )

end


end