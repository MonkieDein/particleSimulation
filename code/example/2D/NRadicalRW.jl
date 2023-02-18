include("../../basics.jl")
# Simulate radical in a particle with random walk
nm = 100
particle = Circle(0,0,nm)

radical_radius = 3
radical_x = particle.min.x + 1e-10 + radical_radius
radical_y = 0

# Multiple Radicals
N = 100
Rad = [Circle(radical_x,radical_y,radical_radius,vx = 100,vy = 0) for r in 1:N]

# video FPS
Δt = (1/60)

# random directional change
Δvt = (1/10)
Δvi = Int(round(Δvt/Δt))

T = collect(0:Δt:30)

video = @animate for (i,t) in ProgressBar(enumerate(T))

    # Plot all radicals 
    plt1 = plot(legend=:outerright , xlabel = "nm",ylabel = "nm",title = "$N instances of radical within Particle",size = (700,600))
    plotShape!(plt1,particle,"particle")
    plotShapes!(plt1,Rad,"")

    L2D = [L2Distance(r.p,particle.p) for r in Rad]
    plt2 = histogram(L2D, label = "Experimental",bins=range(0, nm, length=21),title="Distance from center Distribution",xlims = (0,nm),ylims = (0,N))
    plot(plt1,plt2,layout = (1,2),size = (1400, 600))

    # If time interval equal to radical direction change, then change the radical
    if (i % Δvi)  == 1
        for radical in Rad
            random_direction(radical)
        end
    end

    for radical in Rad
        tempΔt = Δt
        while tempΔt > 0
            colission = checkParticleCollision(particle,radical,tempΔt)
            updateMotion(radical , colission.time * tempΔt)
            if (colission.time < 1)
                radical.v = colission.reflection
            end
            tempΔt -= colission.time * tempΔt
        end
    end
end


gif(video,"animation/NradicalRW.mp4",fps=60)