include("../../reactChem.jl")
include("../../reactStatsPlots.jl")
using StatsBase

wpInitsArray = [0.75] # [0.8,0.95] 
layerRadisArray = [0,60,70] #round.(Int,)
Tgvalues = [25.4] # [25.4,45.1,102.1]
histCounts = []
increments = collect(0:10:20)
N = 1000                                                                                   # N : NumberOfRadical
Vs = fill(0.0,(length(Tgvalues),length(wpInitsArray),N))
histograms = []
for (ntg,tgv) in enumerate(Tgvalues)
for (nwp,wpInit) in enumerate(wpInitsArray)
for incr in increments
# Multi Layer Particle variables
wpEnd = 1.0                                      # initial Wp values
startWpTime , endWpTime = 0.0 , 30*60.0                        # linear Relationship Wp reaction end Time
parRadius = layerRadisArray[length(layerRadisArray)]+incr                           # particle radius
reactionTemp = 70                                               # T : Reaction Temparature (°C)
layerRadius = layerRadisArray[1:(length(layerRadisArray)-1)]    # 2 # (1+nwp)
Tg₀ = fill(tgv,size(layerRadius))                              # initial Tg
Tmon = 106                                                      # monomer temp 106°C

# --- Make particle observable function --- 
layerR = Observable(layerRadius)                                # layerR : layer Radius away from center)
Wp = Node(wpInit)                                               # Wp : Weightage of polymer
layerΔT = Observable(reactionTemp .- Tg₀)                       # layerΔT = T - Tg # (@lift(reactionTemp .- updateTg($Wp,Tg₀,Tmon)) )
layerD = @lift([10^logD($Wp,t,unit="nm") for t in $layerΔT])    # Layers Diffusion coefficient
par = mLparticle(parRadius,Wp,layerR,layerΔT,layerD)                                

# Simulation propagation Time statistics
lPropl = 297                                                                          # Number of propagation steps before simulation end
propTime = @lift(propagationTimeInterval($(par.Wp),reactionTemp))                   # propTime : Time Interval for Monomer to propagate
propStats = getPropStats(wpInit,reactionTemp,endWpTime,lPropl,
                            wpEnd = wpEnd, startTime=startWpTime)                   # get propagation statistics
T = propStats.T                                                                     # T        : Total simulation time interval

# initialize Radicles
zmerInit = 4
maxStepLength = 50 # minimum(diff(par.layerR[])) # 40 #
zmer = Node(zmerInit)
τ =   Node(T/1000) # Node(MinTimeForStepsize(maxStepLength,par.layerD[],wpInit,zmerInit,confident=0.9)  ) #
radRadius = 0.75
Rad = [Radicle([par.obj.radius-radRadius*(1+1e-10),0,0],par,r=radRadius,zmer = zmer,τ = τ) for n in 1:N]
for rad in Rad
    random_position(rad.obj) 
end

# Run simulations
sims = simulate(propTime,Wp,zmer,collect(0:τ[]:T),par,Rad,wpInit,endWpTime,
wpEnd=wpEnd,startTime=startWpTime,wallCollideDo="random",eachStepDo ="random")

simsL2=[L2(pos) for pos in sims.P]
simsLayer=[searchsortedlast(par.layerR[], d ) for d in simsL2]
# Compute the frequency of each value
counts = countmap(reduce(vcat,simsLayer))

# Extract unique values and their frequencies
values = sort(collect(keys(counts)))
frequencies = [counts[value] for value in values]
duration = frequencies ./ sum(frequencies) .* T
barplot = bar(values, duration, xlabel="Layer", ylabel="Duration (second)",ylim=(0,T), label = "", 
title="Histogram of radicles time spend in each Layer",xticks=1:(length(layerRadisArray)-1),xlim=(0.5,length(layerRadisArray)-0.5))
savefig(barplot,wdir("plots/MlayerChgWp/Cliff/LenLimit/")*string(wpInit)*"-wp-Layer-tg-"*string(tgv)*"Radius"*string(parRadius)*".png")

########################################################################
v = floor.(Int,reduce(vcat,simsL2))
num_bins = max(1,Int(floor(maximum(v) - minimum(v))))
counts = countmap(v)
# Extract unique values and their frequencies
values = sort(collect(keys(counts)))
frequencies = [counts[value] for value in values]
duration = frequencies ./ sum(frequencies) .* T
his = bar(values, duration, xlabel="Distance (nm)", ylabel="Duration (second)",ylim=(0,T/10), label = "", 
title="Histogram of radicles time spend in each Distance",xticks=0:5:(layerRadisArray[end]+increments[end]),xlim=(0,(layerRadisArray[end]+increments[end])))
savefig(his,wdir("plots/MlayerChgWp/Cliff/LenLimit/")*string(wpInit)*"-wp-Dist-tg-"*string(tgv)*"Radius"*string(parRadius)*".png")
end
end
end

