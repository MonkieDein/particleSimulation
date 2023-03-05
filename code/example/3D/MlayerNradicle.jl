include("../../statsPlot.jl")

# Multi Layer Particle variables
Wp = 0.95                               # Wp : Weightage of polymer
T = 80                                  # T : Reaction Temparature (°C)
L = DataFrame(                          # L : Layer DataFrame
    R = [0,50,57,62,66],                # R : Radius away from center
    ΔT = (T .- [30,40,60,50,80]) )      # ΔT = T - Tg
par = mLparticle(70,Wp,L)               # particle structure

# Simulation Time and Frame parameter
propTime = propagationTimeInterval(Wp,T)# propTime : Time Interval for Monomer to propagate
lPropl = 8;                             # Number of propagation steps before end
T = lPropl * propTime                   # T        : Total simulation time interval
# τ        : time interval for a Random walk step
maxStepLength = 5
lZmerl = 4;                             # Zmer length (Different for each monomer): 
τ = MinTimeForStepsize(maxStepLength,par.L.D,Wp,lZmerl,confident=0.9)  

# Initialize multiple Radicles variables
N = 100
Rad = [Radicle([-par.obj.radius,0,0],l = length(L.R),zmer = lZmerl,τ = τ) for r in 1:N]
for rad in Rad
    updateRadicle(par,rad)
end

# Monte Carlo Simulation
sims = simulate(τ,propTime,T,par,Rad)

# if you want it to write the simulations result into a folder
# writeSimsData(sims;folder="Data",overwrite=false)

velocityBoxPlot(sims,par) |> display
# savefig("plots/Mlayer3D/velocityBoxPlot.png")

sampleDepthPlot(sims,length(L.R),lPropl,propTime,lZmerl) |> display
# savefig("plots/Mlayer3D/DeepestSamples.png")

