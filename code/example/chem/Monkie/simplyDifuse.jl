include("../../../3Dbasics.jl")
include("../../../chem.jl")
using ProgressBars
using Plots

D = 1e-9 *(1e14)                # D : diffusion constant / coefficient
ℓ = 1/2                         # ℓ : step size  
τ = ℓ/(D)                       # Timestep (τ) : Expected time to take a unit step
T = 1                           # T : Total simulation time interval
Time = collect(0:τ:T)           # Time : Vector of discretized timeframe     
M = length(Time)                # M : Number of steps taken + 1
par = Round(0,0,0,45)           # Define particle initial position and radius
rad = Round(-par.radius,0,0,0)  # Define radicle initial position and radius

# ================================================================= #
# Displacement per τ period x,y,z ∼ N(0,σ²) and E[x²] = 2Dt,
# ⟹ σ² = E[x²] - E[x]² = E[x²] - 0 = E[x²] = 2Dt
# Therefore x,y,z ∼ N(0,σ²) = N(0,2Dt) = √2Dt ⋅ N(0,1) 
# Therefore sqrt(2Dτ)/τ per sec
σxyz = sqrt(2 * D * τ) / τ 
# |v|² = |x|² + |y|² + |z|² = 6Dt ⟹ |v| ∼ √6Dt ⋅ N(0,1) per τ period
# Therefore sqrt(6D/τ) per sec
σv = sqrt(6 * D / τ)
# ================================================================== #

# Identify value to keep track over the time
pos = Vector{coord}(undef, M)   # Radicle position throughout the time

for (i,t) in ProgressBar(enumerate(Time))
    pos[i] = rad.p
    random_velocity(rad,σxyz)
    resampleStepUpdate(par,rad,τ,σxyz)
end

# Calculate distance from center and distance from initial position
locInit = vec(pos[1])
L2center = [L2Distance(vec(par.p),vec(rad_p)) for rad_p in pos]
L2init = [L2Distance(locInit,vec(rad_p)) for rad_p in pos]

# Plot distance from initial position
Plots.plot(Time,L2init)
