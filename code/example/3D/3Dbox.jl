include("../../3Dbasics.jl")
using DynamicalSystems
using GeometryBasics
using GLMakie


box = Square(0,0,0,100)
rad1 = Round(0,0,0,3,vx = 150)

Δt = (1/60)
T = collect(0:Δt:5)

# Create observables so that plot will be update automatically
oBox = Observable( Rect3D(tup(box.min),Tuple(repeat([box.halflen*2],3))) )
oRadical = Observable( Sphere(Point3f0(positionTuple(rad1)),rad1.radius) )
mesh(oBox,color = (:blue,0.3),transparency=true, shading = true)
mesh!(oRadical, color = :red)

for t in ProgressBar(T)
    tempΔt = Δt
    while tempΔt > 0
        colission = checkBoxCollision(box,rad1,tempΔt)
        firstCollisionTime = minimum(vec(colission.time))
        updateMotion(rad1,firstCollisionTime * tempΔt)
        oRadical[] = Sphere(Point3f0(positionTuple(rad1)),rad1.radius)
        if (colission.collide.x + colission.first.x > 1.5)
            rad1.v.x *= -1
        end
        if (colission.collide.y + colission.first.y > 1.5 )
            rad1.v.y *= -1
        end
        if (colission.collide.z + colission.first.z > 1.5 )
            rad1.v.z *= -1
        end
        tempΔt -= firstCollisionTime * tempΔt
    end
    sleep(Δt) # sleep is required for the plot to update in realtime
end





