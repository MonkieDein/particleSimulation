include("../../chem.jl")
using Plots
# Simulate radical in a single direction with fix velocity

function circumference(x,y,radius;n=150)
    t = range(0, 2π, length = n)
    X = x .+ radius .* sin.(t)
    Y = y .+ radius .* cos.(t)
    return(X=X,Y=Y)
end
function plotShape!(plt,shape,label; plotcenter = false)
    plot!(plt,shape.X,shape.Y,label = label)
    if plotcenter
        scatter!(plt,[shape.x],[shape.y],label = label*"-o")
    end
    return plt
end
function plotShapes!(plt,shapes,label; plotcenter = false)
    for shape in shapes
        plotShape!(plt,shape,label, plotcenter = plotcenter)
    end
    return plt
end
# Simulate radical in a particle fixed velocity

# Multi Layer Particle variables
Wp = 0.95    # any random value, zmer is 1, never matter
L = DataFrame(R = [0,20,40,60,80])      # ΔT = T - Tg
par = mLparticle(100,Wp,L,D = (([1,5,0.1,10,0.5]) .^ 2) ./ 2  ) # V²/2 = D ⟹ V = sqrt(2D)

propTime = 10000                          # propTime : Time Interval for Monomer to propagate
τ = 1
T = 1200
Time = collect(0:τ:T)                   # Time     : Vector of discretized timeframe     
N = length(Time)

# Initialize Radicle variables
rad = Radicle([-par.obj.radius,0,0],r=2,l = length(L.R),zmer = 1,τ = τ)
updateRadicle(par,rad)

# create plot shape
particles = [circumference(par.obj.p.x , par.obj.p.y , R) for R in [20,40,60,80,100]]

# position and velocity throughout the time
pos = Vector{coord}(undef, N)      
vel = Vector{coord}(undef, N)     

video = @animate for (i,t) in ProgressBar(enumerate(Time))
    pos[i] = rad.obj.p
    vel[i] = rad.obj.v

    plt = plot(legend=:outerright , xlabel = "nm",ylabel = "nm",title = "Radical within Particle",size = (700,600))
    rad1 = circumference(rad.obj.p.x , rad.obj.p.y , rad.obj.radius)
    plotShapes!(plt,particles,"")
    plotShape!(plt,rad1,"radical")
    multiLbounceStepUpdate(par,rad,τ) 
end
gif(video,"animation/ParticleTest.mp4",fps=60)
# gif(video,"animation/2D/MultiLayer/ParticleTest.mp4",fps=60)

L2center = [L2Distance(vec(par.obj.p),vec(rad_p)) for rad_p in pos]
VoverT = [L2Distance([0.0,0,0],vec(v)) for v in vel]

# Plot distance velocity graph
plt = scatter(VoverT,L2center,xlab="Velocity",ylab="Distance from center")
for r in L.R
    plot!(plt,[minimum(VoverT),maximum(VoverT)],[r,r])
end
plt

# plt2 = plot(Time,L2center,xlab="Time",ylab="Distance from center")
# for r in L.R
#     plot!(plt2,[Time[1],last(Time)],[r,r])
# end
# plt2

#
# ff = findfirst((60 .< L2center .< 80) .& (VoverT .< 1))
# rad.obj.p = pos[ff-1]
# rad.obj.v = vel[ff-1]


