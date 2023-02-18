using ImageMagick
using Setfield
using ProgressBars
using LinearAlgebra
using Distributions
using GeometryBasics
using GLMakie


mutable struct coord
    x :: Float64
    y :: Float64
    z :: Float64
    function coord(x,y,z)
        new(x,y,z)
    end
    function coord(V)
        new(V[1],V[2],V[3])
    end
end

function vec(point::coord)
    return( [point.x,point.y,point.z] )
end

function tup(point::coord)
    return( Tuple(vec(point)) )
end

mutable struct Square
    p :: coord # position of shape
    v :: coord # velocity of shape
    a :: coord # acceleration of shape
    min :: coord # minimum position of the 3D object
    max :: coord # minimum position of the 3D object
    halflen :: Float64 # size of shape
    function Square(x,y,z,halflen;vx=0,vy=0,vz=0,ax=0,ay=0,az=0)
        p = coord(x,y,z)
        new(p,coord(vx,vy,vz),coord(ax,ay,az),coord(vec(p) .- halflen),coord(vec(p) .+ halflen),halflen)
    end
end

function positionTuple(obj)
    return( Tuple(vec(obj.p)) )
end

mutable struct Round 
    p :: coord # position of shape
    v :: coord # velocity of shape
    a :: coord # acceleration of shape
    radius :: Float64 # size of shape
    function Round(x,y,z,radius;vx=0,vy=0,vz=0,ax=0,ay=0,az=0)
        new(coord(x,y,z),coord(vx,vy,vz),coord(ax,ay,az),radius)
    end
end


function DistanceSquare(p1::Vector{Float64},p2::Vector{Float64})
    return sum( (p1 .- p2) .^ 2)
end

function L2Distance(p1::Vector{Float64},p2::Vector{Float64})
    return sqrt( sum( (p1 .- p2) .^ 2) )
end

function random_direction(obj;size = L2Distance(vec(obj.v),zeros(3)))
    ϕ = rand(Uniform(0,2π))
    cosθ = rand(Uniform(-1,1))
    θ = acos( cosθ )
    obj.v.x = size * sin(θ) * cos(ϕ)
    obj.v.y = size * sin(θ) * sin(ϕ)
    obj.v.z = size * cos(θ)
end

function updateMotion(obj,Δt)
    # update position
    obj.p = coord( vec(obj.p) .+ vec(obj.v) .* Δt )
    # update velocity
    obj.v = coord( vec(obj.v) .+ vec(obj.a) .* Δt )
end

function checkBoxCollision(box::Square,ball::Round,Δt)
    
    ball1 = coord( vec(ball.p) .+ vec(ball.v) .* Δt )
    # x = t x₁ + (1 - t) x₀  ; y = t y₁ + (1 - t) y₀ 
    # boxₘᵢₙ.x + ballᵣ = x    ; boxₘᵢₙ.y + ballᵣ = y , substitute x and y and solve for t.
    t = coord(1,1,1)
    if (ball1.x - ball.radius ≤ box.min.x) 
        t.x = min(t.x,(box.min.x + ball.radius - ball.p.x)/(ball1.x - ball.p.x))
    end
    if (ball1.y - ball.radius  ≤ box.min.y) 
        t.y = min(t.y,(box.min.y + ball.radius - ball.p.y)/(ball1.y - ball.p.y))
    end
    if (ball1.z - ball.radius  ≤ box.min.z) 
        t.z = min(t.z,(box.min.z + ball.radius - ball.p.z)/(ball1.z - ball.p.z))
    end
    if (ball1.x + ball.radius  ≥ box.max.x)
        t.x = min(t.x,(box.max.x - ball.radius - ball.p.x)/(ball1.x - ball.p.x))
    end
    if (ball1.y + ball.radius  ≥ box.max.y)
        t.y = min(t.y,(box.max.y - ball.radius - ball.p.y)/(ball1.y - ball.p.y))
    end
    if (ball1.z + ball.radius  ≥ box.max.z) 
        t.z = min(t.z,(box.max.z - ball.radius - ball.p.z)/(ball1.z - ball.p.z))
    end

    return(collide = coord(vec(t) .< 1), first = coord(vec(t) .== minimum(vec(t))), time = t)
end

function checkParticleCollision(particle::Round,ball::Round,Δt)
    ball1 = coord( vec(ball.p) .+ vec(ball.v) .* Δt )

    t = 1
    r = [0,0,0]
    if ( L2Distance(vec(particle.p),vec(ball1)) ≥ (particle.radius - ball.radius) )
        # Solve for t : x = t x₁ + ( 1 - t ) x₀ ;  y = t y₁ + ( 1 - t ) y₀ ;
        # substitute x and y into : (R - r)² = (xₚ - x)² + (yₚ - y)²  and solve for t.   
        a = DistanceSquare(vec(ball.p),vec(ball1)) # (x₀ - x₁)^2 + (y₀ - y₁)^2 + (z₀ - z₁)^2
        b = 2 * sum((vec(particle.p) .- vec(ball.p)) .* (vec(ball.p) .- vec(ball1)))
        c = DistanceSquare(vec(particle.p),vec(ball.p)) - (particle.radius-ball.radius)^2
        if (b^2 - 4*a*c) > 0
            solve = (-b + sqrt(b^2 - 4*a*c))/(2*a)
            if solve < 0
                println(ball)
                error("time < 0, t : ",solve)
            end
            if (0 ≤ solve) && (solve < t)
                t = solve
                # reflection vector r = d - 2(d ⋅ n) n ; 
                # the normal vector n must be normalized ⋅ is the dot product.
                # n assume the center of particle is at position (0,0,0)
                n = normalize( (t .* vec(ball1)) .+ ((1 - t) .* vec(ball.p)) )
                d = vec(ball.v)
                r = d .- 2 * ( dot(n,d) ) .* n
            end 
        end
    end

    return(time = t, reflection = coord(r))
end




