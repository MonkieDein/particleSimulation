include("../../statsPlot.jl")

# Multi Layer Particle variables
Wp = 0.95                               # Wp : Weightage of polymer
T = 70                                  # T : Reaction Temparature (°C)
L = DataFrame(                          # L : Layer DataFrame
    R = [0,50,57,62,66],                # R : Radius away from center
    ΔT = (T .- [30,50,60,70,30]) )      # ΔT = T - Tg
par = mLparticle(70,Wp,L)               # particle structure

# Simulation Time and Frame parameter
propTime = propagationTimeInterval(Wp,T)# propTime : Time Interval for Monomer to propagate
lPropl = 8;                             # Number of propagation steps before end
T = lPropl * propTime                   # T        : Total simulation time interval
# τ        : time interval for a Random walk step
maxStepLength = 10
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

# j = 2 # select a arbitrary monte carlo instance
# # see distancePlot of instance j
# distancePlot(sims,par,j,propTime,lPropl,lZmerl)
# fps = 60
# video = distanceAnim(sims,par,j,propTime,lPropl,lZmerl;videoSec = 20,fps = fps)
# gif(video,"animation/DistTimePlot"*string(j)*".mp4",fps=fps)

# if you want it to write the simulations result into a folder
# writeSimsData(sims;folder="Data",overwrite=false)

velocityBoxPlot(sims,par,τ) |> display
# savefig("plots/Mlayer3D/velocityBoxPlot.png")

sampleDepthPlot(sims,length(L.R),lPropl,propTime,lZmerl) |> display
# savefig("plots/Mlayer3D/DeepestSamples.png")

