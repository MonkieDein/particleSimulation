include("../../basics.jl")
# Simulate radical in a box 

box = Box(0,0,100)
rad1 = Circle(0,0,3,vx = 150,vy = 150)

Δt = (1/60)
T = collect(0:Δt:10)
video = @animate for t in ProgressBar(T)
    plt = plot(legend=:outerright , xlabel = "nm",ylabel = "nm",title = "Radical within Box",size = (700,600))
    plotShape!(plt,box,"box")
    plotShape!(plt,rad1,"radical")
    tempΔt = Δt
    while tempΔt > 0
        colission = checkBoxCollision(box,rad1,tempΔt)
        firstCollisionTime = min(colission.time.x,colission.time.y)
        updateMotion(rad1,firstCollisionTime * tempΔt)
        if (colission.collide.x + colission.first.x > 1.5)
            rad1.v.x *= -1
        end
        if (colission.collide.y + colission.first.y > 1.5 )
            rad1.v.y *= -1
        end
        tempΔt -= firstCollisionTime * tempΔt
    end
end

gif(video,"animation/2D/box.mp4",fps=60)