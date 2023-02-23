include("../../3Dbasics.jl")
include("../../chem.jl")

# Diffusion coefficient parameters
Wp = 1-Wm                               # Wp : Weightage of polymer
T = 80;                                 # T : Reaction Temparature (°C)
Tg_sys = 25;                            # Tg : Glass Transition Temperature (ignore plasticization w/ Monomer)
D = 10^logD(Wp ,T - Tg_sys,unit="nm")   # D : diffusion constant / coefficient
# zmer length and number of propagation 
lZmerl = 4;                             # Zmer length (Different for each monomer): 
lPropl = 14;                            # Number of propagation steps before end

# Simulation Time and Frame parameter
propTime = propagationTimeInterval(Wp,T)   # propTime : Time Interval for Monomer to propagate
fpProp = 300                            # fpProp : frame per propagation
τ = propTime/fpProp                     # τ : time interval for a step
T = lPropl * propTime                   # T : Total simulation time interval
Time = collect(0:τ:T)                   # Time : Vector of discretized timeframe     
N = length(Time)                        # N : Number of steps taken + 1

# Simulation objects
par = Round(0,0,0,50)                   # Define particle initial position and radius
rad = Round(-par.radius,0,0,0)          # Define radicle initial position and radius

# Initialize diffusion Coefficient
Dcurr = D/(lZmerl^(0.5+1.75*Wp))
σxyz = sqrt(2 * Dcurr * τ) / τ  
tnextP = propTime
# Radicle position and propagation Information over the simulation
pos = Vector{coord}(undef, N)           # Radicle position throughout the time
firstProp = DataFrame(Zmer=[lZmerl],D=[Dcurr],t=[0.0])

oParticle = Observable( Sphere(Point3f0(positionTuple(par)),par.radius) )
oRadical = Observable( Sphere(Point3f0(positionTuple(rad)),lZmerl*0.2) )
# GLMakie.wireframe(oParticle)
mesh(oParticle, color = (:blue,0.1),transparency=true, shading = true)
mesh!(oRadical, color = :red)

# Calculate Diffusion coefficient given zmerLength
for (i,t) in ProgressBar(enumerate(Time))
    pos[i] = rad.p
    random_velocity(rad,σxyz)

    tempτ = τ
    while (tnextP - t <= tempτ) # time to propagate
        prePropTime = tnextP - t
        bounceStepUpdate(par,rad,prePropTime) 
        
        tempτ -= prePropTime
        lZmerl += 1
        Dcurr = D/(lZmerl^(0.5+1.75*Wp))
        σxyz = sqrt(2 * Dcurr * τ) / τ  
        push!(firstProp,[lZmerl,Dcurr,t + prePropTime])
        tnextP += propTime
    end
    
    bounceStepUpdate(par,rad,tempτ) 
    
    oRadical[] = Sphere(Point3f0(positionTuple(rad)),lZmerl*0.2)
    sleep(τ) # sleep is required for the plot to update in realtime
end

L2center = [L2Distance(vec(par.p),vec(rad_p)) for rad_p in pos]

maximum(L2center)
minimum(L2center)
