include("../../../3Dbasics.jl")
include("../../../chem.jl")

# Diffusion coefficient parameters
Wm = 0.05                               # Wm : Weightage of monomer
Wp = 1-Wm                               # Wp : Weightage of polymer
T = 80;                                 # T : Reaction Temparature (°C)
Tg_sys = 25;                            # Tg : Glass Transition Temperature (ignore plasticization w/ Monomer)
D = 10^logD(Wp ,T - Tg_sys,unit="nm")   # D : diffusion constant / coefficient
lZmerl = 4;                             # Zmer length (Different for each monomer): 
lPropl = 14;                            # Number of propagation steps before end

# Propagation Time for Monomer
R = 8.314;                              # Gas Constant (J/(K * mol))
Mw = 100                                # Molecule weight (g/mol) ! 
M = (Wm/Mw)*1000;                       # Monomer concentration (mol/g) !
k_p = (2.673e6)*exp(-22.36e3/(R * °C2K(T)));   # Propagation Coefficient ( L/(mol*s) ) !
propTime = 1/(k_p * M);                 # Propagation Time

# Simulation Time and Frame parameter
Dmax = D/((lZmerl)^(0.5+1.75*Wp))       # Maximum Diffusion coefficient
ℓ = 1/2                                 # ℓ : step size  
t₁ = ℓ/(Dmax)                           # t₁ : Minimum time to take a unit step
τ = min( propTime , t₁ )                # τ  : Minimum timestep for a change 
T = lPropl * propTime                   # T  : Total simulation time interval
Time = collect(0:τ:T)                   # Time : Vector of discretized timeframe     
N = length(Time)                        # N : Number of steps taken + 1

# Simulation objects
par = Round(0,0,0,50)                   # Define particle initial position and radius
rad = Round(-par.radius,0,0,0)          # Define radicle initial position and radius

# Initialize diffusion Coefficient
Dcurr = D/(lZmerl^(0.5+1.75*Wp))
tnextP = propTime
# Radicle position and propagation Information over the simulation
pos = Vector{coord}(undef, N)           # Radicle position throughout the time
firstProp = DataFrame(Zmer=[lZmerl],D=[Dcurr],t=[0.0])

# Calculate Diffusion coefficient given zmerLength
for (i,t) in ProgressBar(enumerate(Time))
    if (t >= tnextP) # time to propagate
        lZmerl += 1
        Dcurr = D/(lZmerl^(0.5+1.75*Wp))
        push!(firstProp,[lZmerl,Dcurr,t])
        tnextP += propTime
    end
    pos[i] = rad.p
    σxyz = sqrt(2 * Dcurr * τ) / τ   
    random_velocity(rad,σxyz)
    resampleStepUpdate(par,rad,τ,σxyz)
end

locInit = vec(pos[1])
L2init = [L2Distance(locInit,vec(rad_p)) for rad_p in pos]
L2center = [L2Distance(vec(par.p),vec(rad_p)) for rad_p in pos]

# 问：Discretize Time 的目地
# 问：Wp 会改变吗？
# 问：
