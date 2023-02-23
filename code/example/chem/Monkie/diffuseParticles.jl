include("../../../chem.jl")
using Random
using Plots
using Statistics
using Distributions

# Input Parameter
wₚ = 0.95
Tᵣₓₙ = 70
Tg = 90

# Caculate D
D = 10^logD( wₚ , Tᵣₓₙ - Tg , unit = "nm")

# Simulation Parameter
dt = 0.0000001
NParticles = 100000
NSteps = 1000
Time = collect(0:dt:(NSteps*dt))
particle = zeros(NParticles,3)
Avg_Distance² = zeros(NSteps+1)

# Simulate particle through the time
for (i,t) in enumerate(Time) 
    Avg_Distance²[i] = sum(particle .^2)/NParticles;
    # |x|₂ = sqrt(2 * D * dt) avg axis Displacement per step
    particle = particle + (sqrt(2*D*dt)) .* randn(NParticles,3);
end

# According to slide Avg_Distance² = 6 D t
Plots.plot(Time, Avg_Distance²,label="SampleMean")
Plots.plot!(Time, 6 * D .* Time ,label="ExpectedMean")
