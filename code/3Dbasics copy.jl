using Plots
using ImageMagick
using Setfield
using ProgressBars
using LinearAlgebra
using Distributions

mutable struct coord3D
    x :: Float64
    y :: Float64
    z :: Float64
    function coord3D(x,y,z)
        new(x,y,z)
    end
end

mutable struct Cube
    p :: coord3D # position of shape
    v :: coord3D # velocity of shape
    a :: coord3D # acceleration of shape
    halflen :: Float64 # size of shape
    min :: coord3D # position min
    max :: coord3D # position max
    Surface :: Vector{coord3D} # Surface Vectors
    function Cube(x,y,z,halflen;vx=0,vy=0,vz=0,ax=0,ay=0,az=0)
        xmin = x - halflen
        xmax = x + halflen
        ymin = y - halflen
        ymax = y + halflen
        zmin = z - halflen
        zmax = z + halflen
        X = [ xmin , xmin , xmin , xmin , xmax , xmax , xmax , xmax , xmin ]
        Y = [ ymin , ymin , ymax , ymax , ymax , ymax , ymin , ymin , ymin ]
        Z = [ zmax , zmin , zmin , zmax , zmax , zmin , zmin , zmax , zmax ]
        new(coord3D(x,y,z),coord3D(vx,vy,vz),coord3D(ax,ay,az),halflen,coord3D(xmin,ymin,zmin),coord3D(xmax,ymax,zmax),axisCoord(X,Y,Z))
    end
end

function SphereSurface( x , y , z , radius ; n=26, m=6 )
    cosθ = repeat(range(-1,1, length = m),inner = n)
    ϕ = repeat(range(0, 2π, length = n),outer = m)
    θ = acos.( cosθ )
    X = x .+ radius .* sin.(θ) .* cos.(ϕ)
    Y = y .+ radius .* sin.(θ) .* sin.(ϕ)
    Z = z .+ cos.( θ )
    surface = axisCoord(X,Y,Z)
    return( surface )
end

mutable struct Sphere 
    p :: coord3D # position of shape
    v :: coord3D # velocity of shape
    a :: coord3D # acceleration of shape
    radius :: Float64 # size of shape
    min :: coord3D # position min
    max :: coord3D # position max
    n :: Int # draw discretization
    Surface :: Vector{coord3D} # Y shapes
    function Sphere(x,y,z,radius;vx=0,vy=0,vz=0,ax=0,ay=0,az=0,n=26,m=Int(sqrt(n-1))+1)
        xmin = x - radius
        xmax = x + radius
        ymin = y - radius
        ymax = y + radius
        zmin = z - radius
        zmax = z + radius
        Surface = SphereSurface(x,y,z,radius,n=n,m=m)
        new(coord3D(x,y,z),coord3D(vx,vy,vz),coord3D(ax,ay,az),radius,coord3D(xmin,ymin,zmin),coord3D(xmax,ymax,zmax),n,Surface)
    end
end

box = Cube(0,0,0,100)
Rad = Sphere(0,0,0,3)

plot3d(Xs(box.Surface),Ys(box.Surface),Zs(box.Surface),xlabel ="X",ylabel = "Y",zlabel = "Z")


function L2Distance(p1::coord,p2::coord)
    return sqrt((p1.x - p2.x)^2 + (p1.y - p2.y)^2)
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
                n = normalize([ t * ball1.x + (1 - t) * ball.p.x , t * ball1.y + (1 - t) * ball.p.y ])
                d = [ball.v.x , ball.v.y]
                r = d .- 2 * ( dot(n,d) ) .* n
            end 
        end
        
    end

    return(time = t, reflection = coord(r[1],r[2]))
end




