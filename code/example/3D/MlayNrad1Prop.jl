include("../../statsPlot.jl")

# Multi Layer Particle variables
Wp = 0.95                               # Wp : Weightage of polymer
T = 80                                  # T : Reaction Temparature (°C)
L = DataFrame(                          # L : Layer DataFrame
    R = [0,50,55,60],                # R : Radius away from center
    ΔT = (T .- [25.4,103.3,25.4,103.3]) )      # ΔT = T - Tg
par = mLparticle(65,Wp,L)               # particle structure

# Simulation Time and Frame parameter
propTime = propagationTimeInterval(Wp,T)# propTime : Time Interval for Monomer to propagate
lPropl = 1;                             # Number of propagation steps before end
T = lPropl * propTime                   # T        : Total simulation time interval
# τ        : time interval for a Random walk step
maxStepLength = 10
# P=pos,V=vel,Zmer=zmers,Time=Time,minL=minL
N = 100


P_all = Array{coord}(undef, (0,N))
V_all = Array{coord}(undef, (0,N))
Zmer_all = Array{Int}(undef, (0,N))
minL_all = Array{Int}(undef, (0,N))
Time_all = []
depth_all = []
Zindex = [0]
# Zmer length (Different for each monomer): 
ZmerRange = 4:12
τ = MinTimeForStepsize(maxStepLength,par.L.D,Wp,minimum(ZmerRange),confident=0.9)  
for lZmerl in ZmerRange

    println("Running zmer "*string(lZmerl))

    # Initialize multiple Radicles variables
    Rad = [Radicle([-par.obj.radius,0,0],l = length(L.R),zmer = lZmerl,τ = τ) for r in 1:N]
    for rad in Rad
        updateRadicle(par,rad)
    end

    # Monte Carlo Simulation
    sims = simulate(τ,propTime,T,par,Rad)

    # Saved and combine independent 1 prop simulations
    P_all = vcat(P_all,sims.P)
    V_all = vcat(V_all,sims.V)
    Zmer_all = vcat(Zmer_all,sims.Zmer)
    Time_all = vcat(Time_all,sims.Time)
    minL_all = vcat(minL_all,sims.minL)
    push!(Zindex,length(Time_all))

    ##### J distance plots & Distance animation. J ≤ N
    J = 1
    if J > N
        error("J is bigger than available monte carlo instance N")
    else
        for j in 1:J
            # see distancePlot of instance j
            distancePlot(sims,par,j,propTime,lPropl,lZmerl) |> display
            savefig("plots/Mlayer3D/distPlt_Z"*string(lZmerl)*"_j"*string(j)*".png")
            
            println("Generating Plot Distance Animation ",j)
            video = distanceAnim(sims,par,j,propTime,lPropl,lZmerl;videoSec = 10,fps = 60)
            gif(video,wdir(wdir("animation")*"/Z"*string(lZmerl))*"/DistTimePlotj"*string(j)*".gif",fps=60)
        end
    end
    ###### Write simulation to csv
    # writeSimsData(sims;folder="Data/Z"*string(lZmerl),overwrite=false)
end    
    
Sims_all = (P=P_all,V=V_all,Zmer=Zmer_all,Time=Time_all,minL=minL_all)
    
# velocityBoxPlot(Sims_all,par,τ) |> display
# savefig("plots/Mlayer3D/bindVelocityBoxPlot.png")

depplt = sampleDepthPlot(
    Sims_all,length(L.R),length(ZmerRange),Zindex[2],ZmerRange[1];
    X = 1:length(Sims_all.Time),Z = Zindex  ) |> display
# savefig("plots/Mlayer3D/bindDeepestSamples.png")
