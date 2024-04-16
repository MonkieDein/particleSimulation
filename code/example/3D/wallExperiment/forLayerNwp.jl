include("../../../reactChem.jl")
include("../../../reactStatsPlots.jl")
using StatsBase


wpInitsArray = [0.6916,0.7672,0.8115,0.8408,0.862,0.8797,0.892,0.9034,0.9098]
overallLayerRadisArray = [0 ,64.7, 67.7, 71.85, 75.2, 77.5, 79.7, 81.9, 83.08, 84.72, 87.18 ]
Tgvalues = [25.4,102.1]

result = Dict()
for (nwp,wpInit) in enumerate(wpInitsArray)
    totalLayer =  nwp + 1                                       # respective Layer
    for glossy in 0:totalLayer
        println("running total $totalLayer and glossy at $glossy")
        layerRadisArray = overallLayerRadisArray[1:(totalLayer + 1)]
        TgIndexes = fill(1,totalLayer)
        if glossy != 0 
            TgIndexes[glossy] = 2
        end

        # Multi Layer Particle variables
        colorpalletes = palette(:cool, length(Tgvalues))
        parRadius = layerRadisArray[length(layerRadisArray)]            # particle radius
        reactionTemp = 80                                               # T : Reaction Temparature (°C)
        Tg₀ = Tgvalues[TgIndexes] #  [10,90,20,100,30]                  # initial Tg
        colors = colorpalletes[TgIndexes]
        Tmon = 106                                                      # monomer temp 106°C
        layerRadius = layerRadisArray[1:(length(layerRadisArray)-1)]

        deepest_v = []
        totalWpTime = 30*60.0
        for elapse_min in [15]#0:29
            wpEnd = 1.0                                                     # initial Wp values
            startWpTime , remainWpTime = 0.0 , (totalWpTime - elapse_min*60.0) # linear Relationship Wp reaction end Time
            cur_wp = wpPiecewiseLinear(wpInit,totalWpTime,elapse_min*60.0,wpEnd=wpEnd,startTime=startWpTime)
            N = (30 - elapse_min)*10                                    # N : NumberOfRadical (assume to have linear relationship with remain time)
            # --- Make particle observable function --- 
            layerR = Observable(layerRadius)                                # layerR : layer Radius away from center)
            Wp = Observable(cur_wp)                                               # Wp : Weightage of polymer
            layerΔT = Observable(reactionTemp .- Tg₀)                       # layerΔT = T - Tg # (@lift(reactionTemp .- updateTg($Wp,Tg₀,Tmon)) )
            layerD = @lift([10^logD($Wp,t,unit="nm") for t in $layerΔT])    # Layers Diffusion coefficient
            par = mLparticle(parRadius,Wp,layerR,layerΔT,layerD)     
                                        
            # Simulation propagation Time statistics
            lPropl = 297                                                                          # Number of propagation steps before simulation end
            propTime = @lift(propagationTimeInterval($(par.Wp),reactionTemp))                   # propTime : Time Interval for Monomer to propagate
            propStats = getPropStats(cur_wp,reactionTemp,remainWpTime,lPropl,
                                        wpEnd = wpEnd, startTime=startWpTime)                   # get propagation statistics
            T = propStats.T                                                                     # T        : Total simulation time interval

            # initialize Radicles
            zmerInit = 4
            maxStepLength = 50 # minimum(diff(par.layerR[])) # 40 #
            zmer = Observable(zmerInit)
            τ =   Observable(MinTimeForStepsize(maxStepLength,par.layerD[],cur_wp,zmerInit,confident=0.9)) # Observable(T/10000) #
            radRadius = 0.01
            Rad = [Radicle([par.obj.radius-radRadius*(1+1e-10),0,0],par,r=radRadius,zmer = zmer,τ = τ) for n in 1:N]
            for rad in Rad
                random_position(rad.obj) 
            end
            println("wp $(Wp[]), zmer $(zmer[]), tau $(τ[]), T $T, cur_wp $cur_wp, startTime $startWpTime,remainTime $remainWpTime")
            # Run simulations
            sims = simulate(propTime,Wp,zmer,collect(0:τ[]:T),par,Rad,cur_wp,remainWpTime,
            wpEnd=wpEnd,startTime=startWpTime,wallCollideDo="random",eachStepDo ="random")

            info = getTimeMinMatrix(sims,par.obj.radius,Vars=sims.minDepth)
            v = info.sortedValues[end,:]

            push!(deepest_v,v)
            num_bins = max(1,Int(floor(maximum(v) - minimum(v))))
            his = Plots.histogram(v, normalize=true,bins=num_bins,xticks=0:5:layerRadisArray[end],
            ytick=0:0.1:1.0,ylim=(0,1.0),xlim=(0,overallLayerRadisArray[end]))
            savefig(his,wdir("plots/MlayerChgWp/wallExp/$wpInit-$reactionTemp/total$totalLayer/glossy$glossy")*"/elapse_$elapse_min-DeepestHistogram.png")
        end
        addToDict(result,wpInit,totalLayer,glossy,deepest_v)
        v = reduce(vcat,deepest_v)
        num_bins = max(1,Int(floor(maximum(v) - minimum(v))))
        his = Plots.histogram(v, normalize=true,bins=num_bins,xticks=0:5:layerRadisArray[end],
        ytick=0:0.1:1.0,ylim=(0,1.0),xlim=(0,overallLayerRadisArray[end]))
        savefig(his,wdir("plots/MlayerChgWp/wallExp/$wpInit-$reactionTemp/total$totalLayer/glossy$glossy")*"/overall-DeepestHistogram.png")
    end
end
save_jld("data/wallExperiment.jld",result)

