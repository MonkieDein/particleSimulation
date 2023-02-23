include("../../3Dbasics.jl")
using DynamicalSystems
using GeometryBasics
using GLMakie

particle = Round(0,0,0,109)
rad1 = Round(0,50,0,3,vx = 350)

Δt = (1/60)
T = collect(0:Δt:10)

# Create observables so that plot will be update automatically
oParticle = Observable( Sphere(Point3f0(positionTuple(particle)),particle.radius) )
oRadical = Observable( Sphere(Point3f0(positionTuple(rad1)),rad1.radius) )
mesh(oParticle,color = (:blue,0.3),transparency=true, shading = true)
mesh!(oRadical, color = :red)

for t in ProgressBar(T)
    tempΔt = Δt
    bounceStepUpdate(particle,rad1,tempΔt)
    oRadical[] = Sphere(Point3f0(positionTuple(rad1)),rad1.radius)
    sleep(Δt) # sleep is required for the plot to update in realtime
end
