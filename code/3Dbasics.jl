using ImageMagick
using Setfield
using ProgressBars
using LinearAlgebra
using Distributions
using GeometryBasics


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

function L2(point::coord)
    return( sqrt(sum(vec(point) .^ 2)) )
end

function from3Dto2D(point::coord)
    r = L2(point)
    θ = atan(point.y, point.x)
    return coord( -r .* [cos(θ),sin(θ),0] ) 
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
    function Round(P,radius;V=[0.0,0,0],A=[0.0,0,0])
        new(coord(P),coord(V),coord(A),radius)
    end
    function Round(x,y,z,radius;vx=0,vy=0,vz=0,ax=0,ay=0,az=0)
        new(coord(x,y,z),coord(vx,vy,vz),coord(ax,ay,az),radius)
    end
end


function DistanceSquare(p1::Vector{Float64},p2::Vector{Float64})
    return sum( (p1 .- p2) .^ 2)
end

function L2Distance(p1::Vector{Float64},p2::Vector{Float64})
    return sqrt( DistanceSquare(p1,p2) )
end

function random_position(obj;size = L2Distance(vec(obj.p),zeros(3)))
    ϕ = rand(Uniform(0,2π))
    cosθ = rand(Uniform(-1,1))
    θ = acos( cosθ )
    obj.p =  coord(size .* [ sin(θ)*cos(ϕ) , sin(θ)*sin(ϕ) , cos(θ)])
end

function random_direction(obj;size = L2Distance(vec(obj.v),zeros(3)))
    ϕ = rand(Uniform(0,2π))
    cosθ = rand(Uniform(-1,1))
    θ = acos( cosθ )
    obj.v =  coord(size .* [ sin(θ)*cos(ϕ) , sin(θ)*sin(ϕ) , cos(θ)])
end

function random_velocity(obj,size) # got use
    obj.v = coord(size .* randn(3))
end

function updateMotion(obj,Δt) # got use
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

function collision(particle::Round,ball::Round,Δt)
    ball1 = vec(ball.p) .+ vec(ball.v) .* Δt
    return( L2Distance(vec(particle.p),ball1) ≥ (particle.radius - ball.radius) )
end

# collide 
# : -1 ⟹ from inside.
# : 1 ⟹ from outside.
# : 0 ⟹ center intersection 
# t[1] ⟹ intersection from inner circle
# t[2] ⟹ intersection from outer circle
function intersections(Center, R, ball, Δt;collide=0) # got use
    ball1 =  vec(ball.p) .+ vec(ball.v) .* Δt 
    t = [1.0,1]
    # Solve for t : x = t x₁ + ( 1 - t ) x₀ ;  y = t y₁ + ( 1 - t ) y₀ ;
    # substitute x and y into : (R - r)² = (xₚ - x)² + (yₚ - y)²  and solve for t.   
    # http://paulbourke.net/geometry/circlesphere/
    a = DistanceSquare(vec(ball.p),ball1)
    b = 2 * sum((Center .- vec(ball.p)) .* (vec(ball.p) .- ball1))
    c = DistanceSquare(Center,vec(ball.p)) - (R+(collide * ball.radius))^2
    if (b^2 - 4*a*c) > 0
        solve = (-b + sqrt(b^2 - 4*a*c))/(2*a)
        if (0 ≤ solve) && (solve < t[1]) 
            t[1] = solve
        end
        solve = (-b - sqrt(b^2 - 4*a*c))/(2*a)
        if (0 ≤ solve) && (solve < t[2]) 
            t[2] = solve
        end
    end
    return(t)
end


function closestDistance(P::Vector{Float64}, A::Vector{Float64}, B::Vector{Float64})
    AP = P .- A
    AB = B .- A
    
    dot = sum(AP .* AB)
    lenL = L2Distance(A,B)

    # when AB is just a line with 0 len (i.e. a point), then just return the distance
    if lenL ≤ 0.0
        return L2Distance(A,P)
    end
    t = dot / (lenL^2)
    # since we looking at a line segment we need to clip t to be within A (0) and B (1).
    # C is a point on the line segment that is closest to the point P. 
    C = A .+ AB .* min(max(t,0.0),1.0) 
    return L2Distance(C,P)
end


function collisionNreflection(particle::Round,ball::Round,Δt) # got use
    ball1 =  vec(ball.p) .+ vec(ball.v) .* Δt 
    if ( L2Distance(vec(particle.p),ball1) ≥ (particle.radius - ball.radius) )
        t = intersections(vec(particle.p), particle.radius, ball, Δt,collide = -1)[1]
    else
        t = 1
    end
    r = [0.0,0,0]
    if (t<1)
        # reflection vector r = d - 2(d ⋅ n) n ; 
        # the normal vector n must be normalized ⋅ is the dot product.
        n = normalize( vec(particle.p) .- ((t .* ball1) .+ ((1 - t) .* vec(ball.p))) )
        d = vec(ball.v)
        r = d .- 2 * ( dot(n,d) ) .* n
    end 

    if (t < 0 )
        error("There is error in collision time")
    end
    return(time = t, reflection = coord(r))
end


function bounceStepUpdate(particle::Round,radicle::Round,Δt)  
    tempΔt = Δt
    while tempΔt > 0
        colission = collisionNreflection(particle,radicle,tempΔt)
        updateMotion(radicle,colission.time * tempΔt)
        if (colission.time < 1)
            radicle.v = colission.reflection
        end
        tempΔt -= colission.time * tempΔt
    end
end

function resampleStepUpdate(particle::Round,radicle::Round,Δt,axisSize)
    while (collision(particle,radicle,Δt) == true)
        random_velocity(radicle,axisSize)
    end
    updateMotion(radicle,Δt)
end


