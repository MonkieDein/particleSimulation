include("../../../reactChem.jl")
include("../../../reactStatsPlots.jl")
using StatsBase

wpInitsArray = [0.6916,0.7672,0.8115,0.8408,0.862,0.8797,0.892,0.9034,0.9098]
Tgvalues = [25.4,102.1]
TgIndexes = [1;1;1;1;1;1;1;1;2;1]
colorpalletes = palette(:cool, length(Tgvalues))

# Multi Layer Particle variables
wpInit , wpEnd = 0.8 , 1.0                                      # initial Wp values
startWpTime , endWpTime = 0.0 , 30*60.0                         # linear Relationship Wp reaction end Time
parRadius = layerRadisArray[length(layerRadisArray)]            # 70 # particle radius
reactionTemp = 80                                               # T : Reaction Temparature (°C)
Tg₀ = Tgvalues[TgIndexes] #  [10,90,20,100,30]                  # initial Tg
colors = colorpalletes[TgIndexes]
Tmon = 106                                                      # monomer temp 106°C
layerRadius = layerRadisArray[1:(length(layerRadisArray)-1)]

# --- Make particle observable function --- 
layerR = Observable(layerRadius)                                # layerR : layer Radius away from center)
Wp = Observable(wpInit)                                               # Wp : Weightage of polymer
layerΔT = Observable(reactionTemp .- Tg₀)                       # layerΔT = T - Tg # (@lift(reactionTemp .- updateTg($Wp,Tg₀,Tmon)) )
layerD = @lift([10^logD($Wp,t,unit="nm") for t in $layerΔT])    # Layers Diffusion coefficient
par = mLparticle(parRadius,Wp,layerR,layerΔT,layerD)     
                             
# Simulation propagation Time statistics
lPropl = 5                                                                          # Number of propagation steps before simulation end
propTime = @lift(propagationTimeInterval($(par.Wp),reactionTemp))                   # propTime : Time Interval for Monomer to propagate
propStats = getPropStats(wpInit,reactionTemp,endWpTime,lPropl,
                            wpEnd = wpEnd, startTime=startWpTime)                   # get propagation statistics
T = propStats.T                                                                     # T        : Total simulation time interval

# initialize Radicles
N = 1000                                                                                   # N : NumberOfRadical
zmerInit = 4
maxStepLength = 30 # minimum(diff(par.layerR[])) # 40 #
zmer = Observable(zmerInit)
τ =   Observable(MinTimeForStepsize(maxStepLength,par.layerD[],wpInit,zmerInit,confident=0.9)) # Observable(T/10000) #
radRadius = 0.75
Rad = [Radicle([par.obj.radius-radRadius*(1+1e-10),0,0],par,r=radRadius,zmer = zmer,τ = τ) for n in 1:N]
for rad in Rad
    random_position(rad.obj) 
end

# Run simulations
sims = simulate(propTime,Wp,zmer,collect(0:τ[]:T),par,Rad,wpInit,endWpTime,
wpEnd=wpEnd,startTime=startWpTime,wallCollideDo="random",eachStepDo ="random")

depthHistogramAnim(sims,par,colPalette = :cool,fps=60,videoName=wdir("animation/3D/"*string(wpInit))*"/depthHistogramAnim.mp4",minimum = true,colors=colors)
depthHistogramAnim(sims,par,colPalette = :cool,fps=60,videoName=wdir("animation/3D/"*string(wpInit))*"/positionHistogramAnim.mp4",minimum = false,colors=colors)

anim2D(sims,par,τ,radicalRadius=radRadius,fps=60,videoName=wdir("animation/3D/"*string(wpInit))*"/3DTo2dAnim.mp4",
C = colors)
anim3D(sims,par,radicalRadius=radRadius,fps=60,videoName=wdir("animation/3D/"*string(wpInit))*"/Makie3Danimation.mp4",radC=colors)

dpL = depthPlot(sims,par,length(par.layerR[]),lPropl,propStats,zmerInit;X=sims.Time,Z=[],ntick=2)
savefig(dpL,wdir("plots/MlayerChgWp/"*string(wpInit))*"/depthPlotLayer.png")

dpR = depthPlot(sims,par,0,lPropl,propStats,zmerInit;X=sims.Time,Z=[],ntick=2)
savefig(dpR,wdir("plots/MlayerChgWp/"*string(wpInit))*"/depthPlotRadius.png")

info = getTimeMinMatrix(sims,par.obj.radius,Vars=sims.minDepth)
v = info.sortedValues[size(info.sortedValues)[1],:]
num_bins = max(1,Int(floor(maximum(v) - minimum(v))))
his = Plots.histogram(v, normalize=true,bins=num_bins,xticks=0:5:layerRadisArray[length(layerRadisArray)],
ytick=0:0.1:1.0,ylim=(0,1.0),xlim=(0,layerRadisArray[length(layerRadisArray)]))

savefig(his,wdir("plots/MlayerChgWp/"*string(wpInit))*"/DeepestHistogram.png")



