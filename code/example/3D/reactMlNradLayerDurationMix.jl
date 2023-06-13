include("../../reactChem.jl")
include("../../reactStatsPlots.jl")
using StatsBase

wpInitsArray = [0.8,0.95] 
layerRadisArray = [0,51,57.19,62.11,66.21,69.76,72.89,75.72,78.31,80.68,82.88] #round.(Int,)
Tgvalues = [25.4,45.1,102.1]
colorpalletes = palette(:cool, length(Tgvalues))
TgIndexes = [1;1:3;1:3;1:3]
N = 1000                                                                                   # N : NumberOfRadical
# for (ntg,tgv) in enumerate(Tgvalues)
# ntg = 3
# tgv = Tgvalues[ntg]
for (nwp,wpInit) in enumerate(wpInitsArray)
# nwp = 2
# wpInit = wpInitsArray[nwp]
# Multi Layer Particle variables
wpEnd = 1.0                                      # initial Wp values
startWpTime , endWpTime = 0.0 , 30*60.0                        # linear Relationship Wp reaction end Time
parRadius = layerRadisArray[length(layerRadisArray)]                           # particle radius
reactionTemp = 70                                               # T : Reaction Temparature (°C)
layerRadius = layerRadisArray[1:(length(layerRadisArray)-1)]    # 2 # (1+nwp)
Tg₀ = Tgvalues[TgIndexes] # fill(tgv,size(layerRadius))                              # initial Tg
colors = colorpalletes[TgIndexes]
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
barplot = bar(values, duration, xlabel="Layer", ylabel="Duration (second)",ylim=(0,T), label = "Histogram", 
title="Histogram of radicles time spend in each Layer",xticks=1:(length(layerRadisArray)-1),xlim=(0.5,length(layerRadisArray)-0.5))
for l in 1:length(layerRadius)
    BoxShape(barplot,[l-0.5,l+0.5],[0,T];bw= 1,lc=plot_color(colors[l], 0.1),c = plot_color(colors[l], 0.1))
end
for tgval in 1:length(Tgvalues)
    BoxShape(barplot,[-4,-5],[-4,-5];bw= 1,lc=plot_color(colorpalletes[tgval], 0.1),c = plot_color(colorpalletes[tgval], 0.1),label="Tg: "*string(Tgvalues[tgval])*"°C")
end
savefig(barplot,wdir("plots/MlayerChgWp/Cliff/TimeSpend/")*string(wpInit)*"-wp-Layer-tg-mix.png")
# savefig(barplot,wdir("plots/MlayerChgWp/Cliff/TimeSpend/")*string(wpInit)*"-wp-Layer-tg-"*string(tgv)*".png")

v = floor.(Int,reduce(vcat,simsL2))
num_bins = max(1,Int(floor(maximum(v) - minimum(v))))
counts = countmap(v)

# Extract unique values and their frequencies
values = sort(collect(keys(counts)))
frequencies = [counts[value] for value in values]
duration = frequencies ./ sum(frequencies) .* T
his = bar(values, duration, xlabel="Distance (nm)", ylabel="Duration (second)",ylim=(0,T), label = "Histogram", 
title="Histogram of radicles time spend in each Distance",xticks=0:5:layerRadisArray[length(layerRadisArray)],xlim=(0,layerRadisArray[length(layerRadisArray)]))
for l in 1:length(layerRadius)
    BoxShape(his,layerRadisArray[[l,l+1]],[0,T];bw= 1,lc=plot_color(colors[l], 0.1),c = plot_color(colors[l], 0.1))
end
for tgval in 1:length(Tgvalues)
    BoxShape(his,[-4,-5],[-4,-5];bw= 1,lc=plot_color(colorpalletes[tgval], 0.1),c = plot_color(colorpalletes[tgval], 0.1),label="Tg: "*string(Tgvalues[tgval])*"°C")
end
savefig(his,wdir("plots/MlayerChgWp/Cliff/TimeSpend/")*string(wpInit)*"-wp-Dist-tg-mix.png")


end
# end

