include("../../basics.jl")
# Simulate radical in a particle fixed velocity

particle = Circle(0,0,100)
rad1 = Circle(0,100/sqrt(2),3,vx = 300,vy = 0)

Δt = (1/60)
 
T = collect(0:Δt:10)
video = @animate for t in ProgressBar(T)
    plt = plot(legend=:outerright , xlabel = "nm",ylabel = "nm",title = "Radical within Particle",size = (700,600))
    plotShape!(plt,particle,"particle")
    plotShape!(plt,rad1,"radical")
    tempΔt = Δt
    while tempΔt > 0
        colission = checkParticleCollision(particle,rad1,tempΔt)
        updateMotion(rad1 , colission.time * tempΔt)
        if (colission.time < 1)
            rad1.v = colission.reflection
        end
        tempΔt -= colission.time * tempΔt
    end
end

gif(video,"animation/particle.mp4",fps=60)