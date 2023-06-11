include("../../reactChem.jl")
include("../../reactStatsPlots.jl")

# Multi Layer Particle variables
wpInit , wpEnd = 0.5 , 1.0                                      # initial Wp values
startWpTime , endWpTime = 0.0 , 1.0                            # linear Relationship Wp reaction end Time
parRadius = 70                                                  # particle radius
reactionTemp = 10                                               # T : Reaction Temparature (°C)
Tg₀ = [30,80,50,60,40]                                          # initial Tg
Tmon = 106                                                      # monomer temp 106°C
layerRadius = [0,15,30,50,54]

# --- Make particle observable function --- 
layerR = Observable(layerRadius)                                # layerR : layer Radius away from center)
Wp = Node(wpInit)                                               # Wp : Weightage of polymer
layerΔT = Observable(reactionTemp .- Tg₀)                       # layerΔT = T - Tg # (@lift(reactionTemp .- updateTg($Wp,Tg₀,Tmon)) )
layerD = @lift([10^logD($Wp,t,unit="nm") for t in $layerΔT])    # Layers Diffusion coefficient
par = mLparticle(parRadius,Wp,layerR,layerΔT,layerD)                                

# Simulation propagation Time statistics
lPropl = 499                                                                          # Number of propagation steps before simulation end
propTime = @lift(propagationTimeInterval($(par.Wp),reactionTemp))                   # propTime : Time Interval for Monomer to propagate
propStats = getPropStats(wpInit,reactionTemp,endWpTime,lPropl,
                            wpEnd = wpEnd, startTime=startWpTime)                   # get propagation statistics
T = propStats.T                                                                     # T        : Total simulation time interval

# initialize Radicles
N = 100                                                                                   # N : NumberOfRadical
zmerInit = 2
maxStepLength = 35 # minimum(diff(par.layerR[])) # 40 #
zmer = Node(zmerInit)
radRadius = 0.75
Rad = [Radicle([par.obj.radius-radRadius*1.1,0,0],par,r=radRadius,zmer = zmer,τ = propTime) for n in 1:N]
for rad in Rad
    random_position(rad.obj) 
end

# Run simulations
sims = simulate(propTime,cumsum(propStats.times),par,Rad,wpInit,endWpTime,
wpEnd=wpEnd,startTime=startWpTime,wallCollideDo="random",eachStepDo ="random")

depthHistogramAnim(sims,par,colPalette = :cool,fps=60,videoName="animation/3D/depthHistogramAnimViaProp.mp4")
depthHistogramAnim(sims,par,colPalette = :cool,fps=60,videoName="animation/3D/positionHistogramAnimViaProp.mp4",minimum = false)
