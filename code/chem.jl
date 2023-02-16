using DataFrames
using GLMakie

Region = DataFrame(
    C1 = [ -4.428 ; 26.0 ; 159.0 ; -13.7 ],
    C2 = [ 1.842 ; 37.0 ; 170.0 ; 0.500 ],
    C3 = [ 0 ; 0.0797 ; 0.3664 ; 0 ],
    C4 = [ 8.12e-3 ; 0 ; 0 ; 0 ],
    wpmin = [0 ; 0.795 ; 0.927 ; 0.945],
    wpmax = [0.795 ; 0.927 ; 0.945; 1]
)
# Region is treated as a constant
# ΔT = Tᵣₓₙ - Tg : Tg = 70
function log₁₀D( wₚ , ΔT , Region; scale = 14)
    RegionIndex = searchsortedfirst(Region.wpmax,wₚ)
    C = collect(Region[RegionIndex,["C1","C2","C3","C4"]])
    return( (C[1] - wₚ*C[2]) + (ΔT * C[3]) + (wₚ * ΔT * C[4]) + scale)
end

N = 51
M = 51
Wₚ = collect(range(0,1,N))
ΔT = collect(range(-100,100,M))

logD = [ [log₁₀D( Wₚ[i] ,ΔT[j],Region, scale = 0) for j in 1:M] for i in 1:N ]
logD = reduce(hcat,logD)
lim = FRect3D((-100,0,-16),(200,2,21))

fig = Figure(resolution = (600, 400)) 
ax = Axis3(fig[1,1:2],viewmode = :fit)
# scale!(ax.scene, 1, 1.5, 1)
surface!(ax,ΔT,Wₚ,logD,color = (:blue,0.1),transparency=true)

