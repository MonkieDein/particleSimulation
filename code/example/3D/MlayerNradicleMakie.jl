include("../../chem.jl")
using GLMakie
using Colors

# Multi Layer Particle variables
Wp = 0.95                               # Wp : Weightage of polymer
T = 70                                  # T : Reaction Temparature (°C)
L = DataFrame(                          # L : Layer DataFrame
    R = [0,50,57,62,66],                # R : Radius away from center
    ΔT = (T .- [30,50,60,70,30]) )      # ΔT = T - Tg
par = mLparticle(70,Wp,L)               # particle structure

# Simulation Time and Frame parameter
propTime = propagationTimeInterval(Wp,T)# propTime : Time Interval for Monomer to propagate
lPropl = 4;                             # Number of propagation steps before end
T = lPropl * propTime                   # T        : Total simulation time interval
# τ        : time interval for a Random walk step
maxStepLength = 10
lZmerl = 4;                             # Zmer length (Different for each monomer): 
τ = MinTimeForStepsize(maxStepLength,par.L.D,Wp,lZmerl,confident=0.9)  

# Initialize multiple Radicles variables
N = 100
Rad = [Radicle([-par.obj.radius,0,0],l = length(L.R),zmer = lZmerl,τ = τ) for r in 1:N]
for rad in Rad
    updateRadicle(par,rad)
end

# Monte Carlo Simulation
sims = simulate(τ,propTime,T,par,Rad)
L2CMat = [searchsortedlast(L.R,L2Distance(vec(par.obj.p),vec(rad_p))) for rad_p in sims.P]


radC = range(HSV(-120,1,1), stop=HSV(-360,1,1), length=length(par.L.R)+1)
hexC = [hex(c,:AARRGGBB) for c in radC] # Add alpha 00-FF considering transparency (AA)
C = [parse(Colorant,"0x22"*c[3:8]) for c in hexC ]


oParticle = Sphere(Point3f0(positionTuple(par.obj)),par.obj.radius) 

# GLMakie.wireframe(oParticle)
radplotsize = 1 # sims.Zmer[1,j]*0.4
nRad = N
oRadical = [Observable( Sphere(Point3f0( Tuple(vec(sims.P[1,j])) ),radplotsize) ) for j in 1:nRad]
oRadCol = [Observable( radC[L2CMat[1,j]+1] ) for j in 1:nRad]

fig = mesh(oParticle, color = last(C), shading=true,overdraw=true) #

for l in length(par.L.R):-1:2
    mesh!(Sphere(Point3f0(positionTuple(par.obj)),par.L.R[l]), color = (C[l]),overdraw=true,lightposition = Vec3f0(-150,150,100)) # ,overdraw=true
end

for j in 1:nRad
    mesh!(oRadical[j], color = oRadCol[j] , overdraw=true)
end

fps = 60
secs = 30
I = [1;round.(Int,range(0,1,secs*fps)[Not(1)] .* length(sims.Time))]
frames = 1:length(I)#1:length(sims.Time)
record(fig,"video3D.mp4",frames;framerate = fps) do f #i 
# for (i,t) in ProgressBar(enumerate(sims.Time)) #i in I #
    i = I[f]
    for j in 1:nRad
        oRadCol[j][] = radC[L2CMat[i,j]+1]
        oRadical[j][] = Sphere(Point3f0( Tuple(vec(sims.P[i,j])) ),radplotsize)
    end
    # sleep(1e-10) # sleep is required for the plot to update in realtime
end
