include("../../../chem.jl")
using GLMakie

# Plot Discretization : How fine does each axis has to discretize.
N = 101
M = 101
Wₚ = collect(range(0,1,N)) # collect turn range to a vectors of values
ΔT = collect(range(-100,100,M))

# Calculate logD for each combination of Wₚ and ΔT.
log₁₀D = [ [logD( Wₚ[i] ,ΔT[j],unit = "cm",Region=Region)  for i in 1:N ] for j in 1:M]
log₁₀D = reduce(hcat,log₁₀D) # Turn vector of vector to a 2D Matrices

# Create Figure set resolution, set axis, plot surface.
fig = Figure(resolution = (1200, 800)) 
ax = Axis3(fig[1,1:2],xlabel = "Wₚ",ylabel = "ΔT (°C)",zlabel="log₁₀D (cm²/s)")
surface!(ax,Wₚ,ΔT,log₁₀D,color = (:blue,0.05),transparency=true)
