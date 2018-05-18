# Minmum l2 projection onto the probability simplex. 
# From Fig. 1 in https://stanford.edu/~jduchi/projects/DuchiShSiCh08.pdf
function simplex_projection(v::Vector{Float64}, z::Float64=1.0,
                            minval::Float64=1e-8)
    μ = copy(v);
    sort!(μ, rev=true)
    ρ = maximum(find((μ - (cumsum(μ) - z) ./ collect(1:length(μ))) .> 0))
    θ = (sum(μ[1:ρ]) - z) / ρ
    ret = max.(v - θ, 0)
    ret = max.(ret, minval)
    return ret / sum(ret)
end
