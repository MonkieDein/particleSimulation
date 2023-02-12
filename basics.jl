using Plots
using ImageMagick

struct coordinate
    x :: Float64
    y :: Float64
    function coordinate(x,y)
        new(x,y)
    end
end

struct Box
    p :: coordinate # position of shape
    v :: coordinate # velocity of shape
    a :: coordinate # acceleration of shape
    halflen :: Float64 # size of shape
    xlim :: (Float64,Float64)
    ylim :: (Float64,Float64)
    X :: Vector{Float64} # X shapes
    Y :: Vector{Float64} # Y shapes
    function Box(x,y,halflen;vx=0,vy=0,ax=0,ay=0)
        X = [x - halflen, x - halflen,x + halflen,x + halflen,x - halflen]
        Y = [y - halflen, y + halflen,x + halflen,y - halflen,y - halflen]
        new(coordinate(x,y),coordinate(vx,vy),coordinate(ax,ay),halflen,X,Y)
    end
end

struct Circle 
    p :: coordinate # position of shape
    v :: coordinate # velocity of shape
    a :: coordinate # acceleration of shape
    radius :: Float64 # size of shape

    n :: Int
    X :: Vector{Float64} # X shapes
    Y :: Vector{Float64} # Y shapes
    function Circle(x,y,radius;vx=0,vy=0,ax=0,ay=0,n=150)
        t = range(0, 2π, length = n)
        X = x .+ radius .* sin.(t)
        Y = y .+ radius .* cos.(t)
        new(coordinate(x,y),coordinate(vx,vy),coordinate(ax,ay),radius,n,X,Y)
    end
end

function set_pos(obj,x,y)
    obj.p.x = x
    obj.p.y = y
end

function set_speed(obj,vx,vy)
    obj.v.x = vx
    obj.v.y = vy
end

function set_acceleration(obj,ax,ay)
    obj.a.x = ax
    obj.a.y = ay
end

function updateAll(obj,Δt)
    # update position
    obj.p.x += obj.v.x * Δt
    obj.p.y += obj.v.y * Δt
    # update velocity
    obj.v.x += obj.a.x * Δt
    obj.v.y += obj.a.y * Δt
end

function inCollision(ball::Circle,box::Box)
    
end

function plotShape!(plt,particle,label; plotcenter = false)
    plot!(plt,particle.X,particle.Y,label = label)
    if plotcenter
        scatter!(plt,[particle.x],[particle.y],label = label*"-o")
    end
    return plt
end


particle = Box(0,0,100)
plt = plot(xlabel = "nm",ylabel = "nm",title = "Radical within Particle",size = (600,600))
plotShape!(plt,particle,"particle")

rad1 = Circle(10,50,3,vx = 1,vy = 1)
plotShape!(plt,rad1,"radical")

Δt = (1/100)
T = collect(0:Δt:10)
for t in T

end


