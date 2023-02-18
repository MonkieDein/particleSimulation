include("../../basics.jl")
# Simulate radical in a particle with random walk

particle = Circle(0,0,100)

radical_radius = 3
radical_x = particle.min.x + 1e-10 + radical_radius
radical_y = 0
rad1 = Circle(radical_x,radical_y,radical_radius,vx = 200,vy = 0)

# video FPS
Δt = (1/60)

# random directional change
Δvt = (1/10)
Δvi = Int(round(Δvt/Δt))

T = collect(0:Δt:10)
video = @animate for (i,t) in ProgressBar(enumerate(T))
    plt = plot(legend=:outerright , xlabel = "nm",ylabel = "nm",title = "Radical within Particle",size = (700,600))
    plotShape!(plt,particle,"particle")
    plotShape!(plt,rad1,"radical")
    if (i % Δvi)  == 1
        random_direction(rad1)
    end
    tempΔt = Δt
    while tempΔt > 0
        colission = checkParticleCollision(particle,rad1,tempΔt)
        updateMotion(rad1 , colission.time * tempΔt)
        if (L2Distance(particle.p,rad1.p) > (particle.radius - rad1.radius) )
            println("location error")
        end
        if (colission.time < 1)
            rad1.v = colission.reflection
        end
        tempΔt -= colission.time * tempΔt
    end
end

gif(video,"animation/particleRW.mp4",fps=60)