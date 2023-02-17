using DataFrames
using GLMakie

# Region is treated as a constant
Region = DataFrame(
    C1 = [ -4.428 ; 26.0 ; 159.0 ; -13.7 ],
    C2 = [ 1.842 ; 37.0 ; 170.0 ; 0.500 ],
    C3 = [ 0 ; 0.0797 ; 0.3664 ; 0 ],
    C4 = [ 8.12e-3 ; 0 ; 0 ; 0 ],
)

# ΔT = Tᵣₓₙ - Tg : Tg = 70
# default scale = 14 -> nm²/s , set scale = 0 -> cm²/s .
# logD = ( A - B ⋅ wₚ ) where A = C₁ + ΔT C₃ , and B = C₂ - ΔT C₄
function log₁₀D( wₚ , ΔT , Region; scale = 14)
    A = Region.C1 .+ ΔT .* Region.C3  
    B = Region.C2 .- ΔT .* Region.C4 
    # Intersection of wₚ : (aᵢ - aᵢ₊₁)/(bᵢ - bᵢ₊₁) or (aᵢ₊₁ - aᵢ)/(bᵢ₊₁ - bᵢ)
    RegionMin = [0 ; (A[[1,2,3]] .- A[[2,3,4]]) ./ (B[[1,2,3]] .- B[[2,3,4]])]
    RegionIndex = searchsortedlast(RegionMin,wₚ)
    return( ( A[RegionIndex] - B[RegionIndex] * wₚ ) + scale)
end

# Plot Discretization : How fine does each axis has to discretize.
N = 101
M = 101
Wₚ = collect(range(0,1,N))
ΔT = collect(range(-100,100,M))

# Calculate logD for each combination of Wₚ and ΔT.
logD = [ [log₁₀D( Wₚ[i] ,ΔT[j],Region, scale = 0)  for i in 1:N ] for j in 1:M]
logD = reduce(hcat,logD)

# Create Figure set resolution, set axis, plot surface.
fig = Figure(resolution = (1200, 800)) 
ax = Axis3(fig[1,1:2],xlabel = "Wₚ",ylabel = "ΔT (°C)",zlabel="log₁₀D (cm²/s)")
surface!(ax,Wₚ,ΔT,logD,color = (:blue,0.05),transparency=true)

