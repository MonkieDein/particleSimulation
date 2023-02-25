include("../../chem.jl")
using Colors

# Multi Layer Particle variables
Wp = 0.95                               # Wp : Weightage of polymer
T = 80                                  # T : Reaction Temparature (°C)
L = DataFrame(                          # L : Layer DataFrame
    R = [0,20,40,60,80],                # R : Radius away from center
    ΔT = (T .- [30,80,20,60,30]) )      # ΔT = T - Tg
par = mLparticle(100,Wp,L)              # particle structure

# Simulation Time and Frame parameter
propTime = propagationTimeInterval(Wp,T)# propTime : Time Interval for Monomer to propagate
fpProp = 600                             # fpProp   : frame per propagation
τ = propTime/fpProp                     # τ        : time interval for a Random walk step
lPropl = 14;                            # Number of propagation steps before end
T = lPropl * propTime                   # T        : Total simulation time interval
Time = collect(0:τ:T)                   # Time     : Vector of discretized timeframe     
lTl = length(Time)                        # lTl        : Number of steps taken + 1

# Initialize Radicle variables
lZmerl = 4;                             # Zmer length (Different for each monomer): 
rad = Radicle([-par.obj.radius,0,0],l = length(L.R),zmer = lZmerl,τ = τ)
updateRadicle(par,rad)

# Radicle position and propagation Information over the simulation
pos = Vector{coord}(undef, lTl)           # Radicle position throughout the time
pos = Vector{coord}(undef, lTl)           # Radicle position throughout the time
firstProp = DataFrame(Zmer=[lZmerl],D=[rad.D],t=[0.0])
tnextP = propTime

for (i,t) in ProgressBar(enumerate(Time))
    pos[i] = rad.obj.p
    random_velocity(rad.obj,rad.σx)
    tempτ = τ
    while (tnextP - t <= tempτ) # time to propagate within this time step
        prePropTime = tnextP - t # remaining time before propagation
        multiLbounceStepUpdate(par,rad,prePropTime) 
        
        tempτ -= prePropTime
        rad.zmer += 1
        updateRadicle(par,rad)

        push!(firstProp,[rad.zmer,rad.D,t + prePropTime])
        tnextP += propTime
    end
    
    multiLbounceStepUpdate(par,rad,tempτ) 
    
end

oParticle = Sphere(Point3f0(positionTuple(par.obj)),par.obj.radius) 
oRadical = Observable( Sphere(Point3f0(positionTuple(rad.obj)),rad.zmer*0.4) )
# GLMakie.wireframe(oParticle)
mesh(oParticle, color = (:grey), shading=true,overdraw=true,lightposition = Vec3f0(-150,150,100)) #
C = range(HSV(0,1,1), stop=HSV(-360,1,1), length=length(par.L.R)+1);
for l in length(par.L.R):-1:2
    mesh!(Sphere(Point3f0(positionTuple(par.obj)),par.L.R[l]), color = (C[l]),overdraw=true,lightposition = Vec3f0(-150,150,100)) # ,overdraw=true
end
mesh!(oRadical, color = last(C) , overdraw=true)

oRadical[] = Sphere(Point3f0(positionTuple(rad.obj)),rad.zmer*0.4)
sleep(τ) # sleep is required for the plot to update in realtime
