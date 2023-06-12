include("../../reactChem.jl")
include("../../reactStatsPlots.jl")
using StatsBase

wpInitsArray = [0.783,0.822,0.858,0.887,0.893,0.909,0.919,0.926,0.936] # LinRange(0.70,0.99,4)
layerRadisArray = [0,51,57.19,62.11,66.21,69.76,72.89,75.72,78.31,80.68,82.88] #round.(Int,)
Tgvalues = [25.4,45.1,102.1]
histCounts = []
N = 1000                                                                                   # N : NumberOfRadical
Vs = fill(0.0,(length(Tgvalues),length(wpInitsArray),N))
histograms = []
for (ntg,tgv) in enumerate(Tgvalues)
for (nwp,wpInit) in enumerate(wpInitsArray)
    # Multi Layer Particle variables
    wpEnd = 1.0                                      # initial Wp values
    startWpTime , endWpTime = 0.0 , 30*60.0                        # linear Relationship Wp reaction end Time
    parRadius = layerRadisArray[(2+nwp)]                           # particle radius
    reactionTemp = 70                                               # T : Reaction Temparature (°C)
    layerRadius = layerRadisArray[1:2] # (1+nwp)
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
    
    info = getTimeMinMatrix(sims,par.obj.radius,Vars=sims.minDepth)

    v = info.sortedValues[size(info.sortedValues)[1],:]
    num_bins = max(1,Int(floor(maximum(v) - minimum(v))))
    his = Plots.histogram(v, normalize=true,bins=num_bins,xticks=0:5:layerRadisArray[length(layerRadisArray)],
    ytick=0:0.1:1.0,ylim=(0,1.0),xlim=(0,layerRadisArray[length(layerRadisArray)]), xlabel="Distance (nm) from center", 
    ylabel="Probability", title="Histogram of closest distance from center reach by radicles", legend=false)
    Vs[ntg,nwp,:] = info.sortedValues[size(info.sortedValues)[1],:]

    savefig(his,wdir("plots/MlayerChgWp/Cliff/DeepestHistogram/")*string(wpInit)*"-wp-tg-"*string(tgv)*".png")

end
end

using DelimitedFiles
writedlm( "data/"*string(Tgvalues[1])*".csv",  Vs[1,:,:], ',')
writedlm( "data/"*string(Tgvalues[2])*".csv",  Vs[2,:,:], ',')
writedlm( "data/"*string(Tgvalues[3])*".csv",  Vs[3,:,:], ',')

# v = Vs[3,9,:]
# num_bins = Int(floor(maximum(v) - minimum(v)))
# Plots.histogram(v, normalize=true,bins=num_bins,xticks=0:5:layerRadisArray[length(layerRadisArray)],
# ytick=0:0.1:1.0,ylim=(0,1.0),xlim=(0,layerRadisArray[length(layerRadisArray)]))
