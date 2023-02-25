include("../../chem.jl")
using Plots

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
par = mLparticle(100,Wp,L,D = (([3,15,0.5,30,1.5]) .^ 2) ./ 2  ) # V²/2 = D ⟹ V = sqrt(2D)

propTime = 10000                          # propTime : Time Interval for Monomer to propagate
τ = 1
T = 1200
Time = collect(0:τ:T)                   # Time     : Vector of discretized timeframe     
lTl = length(Time)

# Initialize Radicle variables
N = 100
radr = 2
Rad = [Radicle([-par.obj.radius+radr,0,0],r=radr,l = length(L.R),zmer = 1,τ = τ) for r in 1:N]
for rad in Rad
    updateRadicle(par,rad)
end
# create plot shape
particles = [circumference(par.obj.p.x , par.obj.p.y , R) for R in [L.R[Not(1)];par.obj.radius]]

# position and velocity throughout the time
# pos = Array{coord}(undef, (lTl,N))      
# vel = Array{coord}(undef, (lTl,N))     
# levels = zeros(Int, (lTl,N))     

video = @animate for (i,t) in ProgressBar(enumerate(Time)) # enumerate(Time) #
    plt = plot(legend=:outerright , xlabel = "nm",ylabel = "nm",title = "Radical within multi Layer Particle",size = (700,600))
    radicles = [circumference(rad.obj.p.x , rad.obj.p.y , rad.obj.radius) for rad in Rad]
    plotShapes!(plt,particles,"")
    plotShapes!(plt,radicles,"")
    for rad in Rad
        size = L2Distance(vec(rad.obj.v),[0.0,0,0])
        dir = rand(Uniform(0,2π))
        rad.obj.v.x = size .* sin.(dir)
        rad.obj.v.y = size .* cos.(dir)
    end
    # pos[i,:] = [rad.obj.p for rad in Rad]
    # vel[i,:] = [rad.obj.v for rad in Rad]
    # levels[i,:] = [rad.l for rad in Rad]
    for rad in Rad
        multiLbounceStepUpdate(par,rad,τ) 
    end
end

gif(video,"animation/2D/MultiLayer/Nradicle.mp4",fps=60)

L2CMat = [L2Distance(vec(par.obj.p),vec(rad_p)) for rad_p in pos]
VtMat  = [L2Distance([0.0,0,0],vec(v)) for v in vel]

L2center = reduce(vcat,L2CMat)
VoverT = reduce(vcat,VtMat)

# Plot distance velocity graph
plt = scatter(VoverT,L2center,xlab="Velocity",ylab="Distance from center")
for r in L.R
    plot!(plt,[minimum(VoverT),maximum(VoverT)],[r,r])
end
plt


# ff = findfirst( (80 .< L2CMat .< 100) .& (VtMat .> 10) )

# in_P = vec(pos[ff[1]-1,ff[2]])
# in_V = vec(vel[ff[1]-1,ff[2]])
# in_L = levels[ff[1]-1,ff[2]]

# radicle = Radicle(in_P,V = in_V,r=radr,l = in_L,zmer = 1,τ = τ,σx = L2(coord(in_V)) )
# updateRadicle(par,radicle)
# radicle.obj.p
# radicle.obj.v

# L2(pos[ff[1]-1,ff[2]])
# L2(vel[ff[1]-1,ff[2]])

# pos[ff]
# vel[ff]

# L2(pos[ff])
# L2(vel[ff])
