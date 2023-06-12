include("../../reactChem.jl")
include("../../reactStatsPlots.jl")
using StatsBase

# Multi Layer Particle variables
wpInit , wpEnd = 0.8 , 1.0                                      # initial Wp values
startWpTime , endWpTime = 0.0 , 0.02                            # linear Relationship Wp reaction end Time
parRadius = 70                                                  # particle radius
reactionTemp = 70                                               # T : Reaction Temparature (°C)
Tg₀ = [10,90,20,100,30]                                          # initial Tg
Tmon = 106                                                      # monomer temp 106°C
layerRadius = [0,15,30,50,54]

# --- Make particle observable function --- 
layerR = Observable(layerRadius)                                # layerR : layer Radius away from center)
Wp = Node(wpInit)                                               # Wp : Weightage of polymer
layerΔT = Observable(reactionTemp .- Tg₀)                       # layerΔT = T - Tg # (@lift(reactionTemp .- updateTg($Wp,Tg₀,Tmon)) )
layerD = @lift([10^logD($Wp,t,unit="nm") for t in $layerΔT])    # Layers Diffusion coefficient
par = mLparticle(parRadius,Wp,layerR,layerΔT,layerD)     
                             
# Simulation propagation Time statistics
lPropl = 6                                                                          # Number of propagation steps before simulation end
propTime = @lift(propagationTimeInterval($(par.Wp),reactionTemp))                   # propTime : Time Interval for Monomer to propagate
propStats = getPropStats(wpInit,reactionTemp,endWpTime,lPropl,
                            wpEnd = wpEnd, startTime=startWpTime)                   # get propagation statistics
T = propStats.T                                                                     # T        : Total simulation time interval

# initialize Radicles
N = 100                                                                                   # N : NumberOfRadical
zmerInit = 2
maxStepLength = 35 # minimum(diff(par.layerR[])) # 40 #
zmer = Node(zmerInit)
τ =   Node(T/12000) # Node(MinTimeForStepsize(maxStepLength,par.layerD[],wpInit,zmerInit,confident=0.9)) # Node(T/10000) #
radRadius = 0.75
Rad = [Radicle([par.obj.radius-radRadius*(1+1e-10),0,0],par,r=radRadius,zmer = zmer,τ = τ) for n in 1:N]
for rad in Rad
    random_position(rad.obj) 
end

# Run simulations
sims = simulate(propTime,Wp,zmer,collect(0:τ[]:T),par,Rad,wpInit,endWpTime,
wpEnd=wpEnd,startTime=startWpTime,wallCollideDo="random",eachStepDo ="random")

depthHistogramAnim(sims,par,colPalette = :cool,fps=60,videoName=wdir("animation/3D/"*string(wpInit))*"/depthHistogramAnim.mp4")
depthHistogramAnim(sims,par,colPalette = :cool,fps=60,videoName=wdir("animation/3D/"*string(wpInit))*"/positionHistogramAnim.mp4",minimum = false)

anim2D(sims,par,τ,radicalRadius=radRadius,fps=60,videoName=wdir("animation/3D/"*string(wpInit))*"/3DTo2dAnim.mp4")
anim3D(sims,par,radicalRadius=radRadius,fps=60,videoName=wdir("animation/3D/"*string(wpInit))*"/Makie3Danimation.mp4")

dpL = depthPlot(sims,par,length(par.layerR[]),lPropl,propStats,zmerInit;X=sims.Time,Z=[],ntick=2)
savefig(dpL,wdir("plots/MlayerChgWp/"*string(wpInit))*"/depthPlotLayer.png")

dpR = depthPlot(sims,par,0,lPropl,propStats,zmerInit;X=sims.Time,Z=[],ntick=2)
savefig(dpR,wdir("plots/MlayerChgWp/"*string(wpInit))*"/depthPlotRadius.png")

velBoxPlot = velocityBoxPlot(sims,par,propStats,reactionTemp,Tg₀,τ,
    center=[0.0,0,0],bar_width = 0.8,qs = [0.05,0.25,0.75,0.95],Tmon=106)
savefig(velBoxPlot,wdir("plots/MlayerChgWp/"*string(wpInit))*"/velocityBoxPlot.png")

plts = velocityScatter(sims)
savefig(plts[1],wdir("plots/MlayerChgWp/"*string(wpInit))*"/velocityScatterProp1.png")

info = getTimeMinMatrix(sims,par.obj.radius,Vars=sims.minDepth)
v = info.sortedValues[size(info.sortedValues)[1],:]
num_bins = max(1,Int(floor(maximum(v) - minimum(v))))
his = Plots.histogram(v, normalize=true,bins=num_bins,xticks=0:5:layerRadisArray[length(layerRadisArray)],
ytick=0:0.1:1.0,ylim=(0,1.0),xlim=(0,layerRadisArray[length(layerRadisArray)]))

savefig(his,wdir("plots/MlayerChgWp/"*string(wpInit))*"/DeepestHistogram.png")



