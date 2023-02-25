using Plots
using ImageMagick
using Setfield
using ProgressBars
using LinearAlgebra
using Distributions

mutable struct coord
    x :: Float64
    y :: Float64
    function coord(x,y)
        new(x,y)
    end
end

mutable struct Box
    p :: coord # position of shape
    v :: coord # velocity of shape
    a :: coord # acceleration of shape
    halflen :: Float64 # size of shape
    min :: coord # position min
    max :: coord # position max
    X :: Vector{Float64} # X shapes
    Y :: Vector{Float64} # Y shapes
    function Box(x,y,halflen;vx=0,vy=0,ax=0,ay=0)
        xmin = x - halflen
        xmax = x + halflen
        ymin = y - halflen
        ymax = y + halflen
        X = [ xmin , xmin , xmax , xmax , xmin ]
        Y = [ ymin , ymax , ymax , ymin , ymin ]
        new(coord(x,y),coord(vx,vy),coord(ax,ay),halflen,coord(xmin,ymin),coord(xmax,ymax),X,Y)
    end
end

mutable struct Circle 
    p :: coord # position of shape
    v :: coord # velocity of shape
    a :: coord # acceleration of shape
    radius :: Float64 # size of shape
    min :: coord # position min
    max :: coord # position max
    n :: Int # draw discretization
    X :: Vector{Float64} # X shapes
    Y :: Vector{Float64} # Y shapes
    function Circle(x,y,radius;vx=0,vy=0,ax=0,ay=0,n=150)
        xmin = x - radius
        xmax = x + radius
        ymin = y - radius
        ymax = y + radius
        arc = circumference(x,y,radius,n=n)
        new(coord(x,y),coord(vx,vy),coord(ax,ay),radius,coord(xmin,ymin),coord(xmax,ymax),n,arc.X,arc.Y)
    end
end

function L2Distance(p1::coord,p2::coord)
    return sqrt((p1.x - p2.x)^2 + (p1.y - p2.y)^2)
end

function circumference(x,y,radius;n=150)
    t = range(0, 2π, length = n)
    X = x .+ radius .* sin.(t)
    Y = y .+ radius .* cos.(t)
    return(X=X,Y=Y)
end

function set_pos(obj,x,y)
    obj.p.x = x
    obj.p.y = y
end

function set_speed(obj,vx,vy)
    obj.v.x = vx
    obj.v.y = vy
end

function random_direction(obj)
    size = L2Distance(obj.v,coord(0,0))
    dir = rand(Uniform(0,2π))
    obj.v.x = size .* sin.(dir)
    obj.v.y = size .* cos.(dir)
end

function set_acceleration(obj,ax,ay)
    obj.a.x = ax
    obj.a.y = ay
end

function updateMotion(obj,Δt)
    # update position
    obj.p.x += obj.v.x * Δt
    obj.p.y += obj.v.y * Δt
    # update limits
    obj.min.x = obj.p.x - obj.radius
    obj.max.x = obj.p.x + obj.radius
    obj.min.y = obj.p.y - obj.radius
    obj.max.y = obj.p.y + obj.radius
    # update velocity
    obj.v.x += obj.a.x * Δt
    obj.v.y += obj.a.y * Δt
    # update perimeter location
    arc = circumference( obj.p.x , obj.p.y , obj.radius , n=obj.n )
    obj.X = arc.X
    obj.Y = arc.Y
end

function checkBoxCollision(box::Box,ball::Circle,Δt)
    ball1 = coord(ball.p.x + ball.v.x * Δt , ball.p.y + ball.v.y * Δt)
    # x = t x₁ + (1 - t) x₀  ; y = t y₁ + (1 - t) y₀ 
    # boxₘᵢₙ.x + ballᵣ = x    ; boxₘᵢₙ.y + ballᵣ = y , substitute x and y and solve for t.
    t = coord(1,1)
    if (ball1.x - ball.radius ≤ box.min.x) 
        t.x = min(t.x,(box.min.x + ball.radius - ball.p.x)/(ball1.x - ball.p.x))
    end
    if (ball1.x + ball.radius  ≥ box.max.x)
        t.x = min(t.x,(box.max.x - ball.radius - ball.p.x)/(ball1.x - ball.p.x))
    end
    if (ball1.y - ball.radius  ≤ box.min.y) 
        t.y = min(t.y,(box.min.y + ball.radius - ball.p.y)/(ball1.y - ball.p.y))
    end
    if (ball1.y + ball.radius  ≥ box.max.y)
        t.y = min(t.y,(box.max.y - ball.radius - ball.p.y)/(ball1.y - ball.p.y))
    end

    return(collide = coord(t.x < 1,t.y < 1), first = coord(t.x ≤ t.y,t.y ≤ t.x), time = t)
end

function plotShapes!(plt,shapes,label; plotcenter = false)
    for shape in shapes
        plotShape!(plt,shape,label, plotcenter = plotcenter)
    end
    return plt
end


function plotShape!(plt,shape,label; plotcenter = false)
    plot!(plt,shape.X,shape.Y,label = label)
    if plotcenter
        scatter!(plt,[shape.x],[shape.y],label = label*"-o")
    end
    return plt
end


function checkParticleCollision(particle::Circle,ball::Circle,Δt)
    ball1 = coord(ball.p.x + ball.v.x * Δt , ball.p.y + ball.v.y * Δt)

    t = 1
    r = [0,0]
    if ( L2Distance(particle.p,ball1) ≥ (particle.radius - ball.radius) )
        # Solve for t : x = t x₁ + ( 1 - t ) x₀ ;  y = t y₁ + ( 1 - t ) y₀ ;
        # substitute x and y into : (R - r)² = (xₚ - x)² + (yₚ - y)²  and solve for t.   
        a = (ball.p.x - ball1.x)^2 + (ball.p.y - ball1.y)^2
        b = 2*((particle.p.x - ball.p.x)*(ball.p.x - ball1.x) + (particle.p.y - ball.p.y)*(ball.p.y - ball1.y))
        c = (particle.p.x - ball.p.x)^2 + (particle.p.y - ball.p.y)^2 - (particle.radius-ball.radius)^2
        if (b^2 - 4*a*c) > 0
            solve = (-b + sqrt(b^2 - 4*a*c))/(2*a)
            if solve < 0
                println("time error < 0, t : ",solve)
                println(ball)
            end
            if (0 ≤ solve) && (solve < t)
                t = solve
                # reflection vector r = d - 2(d ⋅ n) n ; 
                # the normal vector n must be normalized ⋅ is the dot product.
                n = normalize([particle.p.x - (t * ball1.x + (1 - t) * ball.p.x) ,particle.p.y -  (t * ball1.y + (1 - t) * ball.p.y) ])
                d = [ball.v.x , ball.v.y]
                r = d .- 2 * ( dot(n,d) ) .* n
            end 
        end
        
    end

    return(time = t, reflection = coord(r[1],r[2]))
end




