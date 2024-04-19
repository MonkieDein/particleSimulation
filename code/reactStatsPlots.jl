include("reactChem.jl")
using Plots
using Colors
using CSV
using Printf

function circumference(x,y,radius;n=150)
    t = range(0, 2π, length = n)
    X = x .+ radius .* sin.(t)
    Y = y .+ radius .* cos.(t)
    return(X=X,Y=Y)
end

function plotShapes!(plt,shapes,label; plotcenter = false,colors = nothing,lcs=colors)
    if isnothing(colors)
        for shape in shapes
            plotShape!(plt,shape,label, plotcenter = plotcenter)
        end
    else
        for (shape,color,lc) in zip(shapes,colors,lcs)
            plotShape!(plt,shape,label, plotcenter = plotcenter,color=color,lc=lc)
        end
    end
    return plt
end


function plotShape!(plt,shape,label; plotcenter = false,color = nothing,lc=color)
    if isnothing(color)
        Plots.plot!(plt,shape.X,shape.Y,label = label)
    else
        Plots.plot!(plt,shape.X,shape.Y,seriestype=[:shape,],c=color,lc=lc,label = label)
    end
    if plotcenter
        Plots.scatter!(plt,[shape.x],[shape.y],label = label*"-o")
    end
    return plt
end

function anim2D(sims,par,τ;radicalRadius = 0.5,fps=60,secs=length(sims.Time)/fps,
    videoName="animation/3D/3DTo2dAnim.mp4",
    radC = palette(:linear_kryw_0_100_c71_n256, size(sims.P)[2]), 
    C = [parse(Colorant,"0x66"*c[3:8]) for c in [hex(c,:AARRGGBB) for c in palette(:cool, par.layerR[])] ])
    
    N = size(sims.P)[2]
    pos2D = [from3Dto2D(p) for p in sims.P]
    particles = [circumference(par.obj.p.x , par.obj.p.y , R) for R in reverse([par.layerR[][Not(1)];par.obj.radius])]
    I =  round.(Int,range(0,1,round(Int,secs*fps)+1)[Not(1)] .* length(sims.Time))
    video = @animate for i in ProgressBar(I) # 
        plt = Plots.plot(legend=:outerright , xlim=(-par.obj.radius,par.obj.radius), ylim=(-par.obj.radius,par.obj.radius), 
        xlabel = "nm",ylabel = "nm",title = "Radical within multi Layer Particle",size = (700,600))
        Plots.annotate!(plt,[ (par.obj.radius/1.1,par.obj.radius/1.01,(@sprintf("Time: %.7f sec",sims.Time[i]),10,:black)) ] )
        Plots.annotate!(plt,[ (par.obj.radius/1.1,par.obj.radius/1.07,("Oligomer: "*string(sims.Zmer[i,1]),10,:black)) ] )
        radicles = [circumference(pos2D[i,j].x , pos2D[i,j].y , radicalRadius, n=16) for j in 1:N]
        plotShapes!(plt,particles,"",colors=reverse(C))
        plotShapes!(plt,radicles,"",colors=radC)
    end
    mp4(video,videoName,fps=fps)
    return nothing
end

function anim3D(sims,par;radicalRadius = 0.5,fps=60,secs=length(sims.Time)/fps,
    videoName="animation/3D/Makie3Danimation.mp4",
    radC = palette(:cool, par.layerR[]), #range(HSV(-120,1,1), stop=HSV(-360,1,1), length=length(par.layerR[]))
    C = [parse(Colorant,"0x33"*c[3:8]) for c in [hex(c,:AARRGGBB) for c in radC] ])
    
    N = size(sims.P)[2]
    L2CMat = [searchsortedlast(par.layerR[],L2Distance(vec(par.obj.p),vec(rad_p))) for rad_p in sims.P]
    # make 3D objects
    oParticle = Sphere(Point3f0(positionTuple(par.obj)),par.obj.radius) 
    oRadical = [Observable( Sphere(Point3f0( Tuple(vec(sims.P[1,j])) ),radicalRadius) ) for j in 1:N]
    oRadCol = [Observable( radC[L2CMat[1,j]] ) for j in 1:N]


    # Initialize first frame, then use observable to update
    lightPos = Vec3f0( 50,-180,10 )

    fig = mesh(oParticle, color = last(C), shading=true,
    overdraw=true,lightposition = lightPos ,
    show_axis=false,figure = (resolution = (1600, 1600),)) #
    # ,transparency = true
    # display(fig)
    for l in length(par.layerR[]):-1:2
        mesh!(Sphere(Point3f0(positionTuple(par.obj)),par.layerR[][l]), color = (C[l-1]),overdraw=true,lightposition = lightPos ) 
    end
    for j in 1:N
        mesh!(oRadical[j], color = oRadCol[j] , overdraw=true,lightposition = lightPos)
    end

    I =  round.(Int,range(0,1,round(Int,secs*fps)+1)[Not(1)] .* length(sims.Time))
    record(fig,videoName,I;framerate = fps) do i
    # for i in ProgressBar(I)
        for j in 1:N
            oRadCol[j][] = radC[L2CMat[i,j]]
            oRadical[j][] = Sphere(Point3f0( Tuple(vec(sims.P[i,j])) ),radicalRadius)
        end
        # sleep(1e-10) # sleep is required if plot fail to update in realtime
    end
    return nothing
end

function getTimeMinMatrix(sims,maxVars;Z=[],Vars=sims.minL)
    lTl = length(sims.Time)
    N = size(sims.P,2)
    values = ones(typeof(Vars[1]),(lTl,N)) .+ maxVars
    sortedValues = ones(typeof(Vars[1]),(lTl,N)) .+ maxVars
    if isempty(Z)
        for i in 1:lTl
            values[i,:] = min.(values[max(i-1,1),:],Vars[i,:])
            sortedValues[i,:] = sort(copy(values[i,:]))
        end
    else
        n = length(Z)
        for z_i in 2:n
            for i in (Z[z_i-1]+1):(Z[z_i])
                values[i,:] = min.(values[max(i-1,Z[z_i-1]+1),:],Vars[i,:])
                sortedValues[i,:] = sort(copy(values[i,:]))
            end
        end
    end
    return (values=values,sortedValues=sortedValues)
end


function depthPlot(sims,par,lLl,lPropl,propStats,lZmerl;X=sims.Time,Z=[],ntick=1,
    colors = palette(:cool, lLl))
    N = size(sims.P,2)
    propTimes = cumsum(propStats.times) 

    if lLl == 0 # Radius depth plot
        colors = palette(:cool)
        minDepthInfo = getTimeMinMatrix(sims,par.obj.radius,Vars=sims.minDepth)
        pltProb = Plots.heatmap(X, 1:N,transpose(minDepthInfo.sortedValues),color = colors,colorbar_title="radius",
            title="Distribution of Deepest Radius reached",legend=:outerright, cbar=:outerright,xlabel="Time",
            ylabel="No of samples",ylims=(0,N*1.1) , dpi=300, clims=(0, par.obj.radius)
        )
    else # Layer depth plots
        layerInfo = getTimeMinMatrix(sims,lLl,Z=Z)
        existLayers = sort(unique(layerInfo.sortedValues))
        # Plot Heat Map
        pltProb = Plots.heatmap(X, 1:N,transpose(layerInfo.sortedValues),
            color = ( length(existLayers) > 1 ? colors[existLayers] : repeat(colors[existLayers],2)),
            title="Distribution of Deepest Layer reached",legend=:outerright,xlabel="Time",
            ylabel="No of samples",ylims=(0,N*1.1),colorbar=false,leg_title ="Layer" , dpi=300
        )
        # Plot Label
        for (i,c) in enumerate(colors)
            Plots.scatter!(pltProb,[],[],color=c,label=string(i))
        end
    end
    # Add customize axis label
    if !isempty(Z)
        Plots.plot!(pltProb,xticks = [])
        zn = length(Z)
        for z_i in 2:zn
            XR = round.(Int,range(Z[z_i-1]+1,Z[z_i],ntick+2))
            for x in XR[2:2+ntick-1]
                Plots.plot!(pltProb,[x,x],[0,1],color=:black,label="")
                Plots.annotate!(pltProb,[ (x,N*(-.02),(string(round(sims.Time[x],sigdigits=2)),5,:black)) ] )
            end
        end
    end
    # Add Zmer Length
    for p in 1:lPropl
        t_i = propTimes[p]
        Plots.plot!(pltProb,[t_i,t_i],[0,N],lc=:black,label="")
        Plots.annotate!(pltProb,[(t_i-0.5*propStats.times[p],N*1.05,(string(lZmerl+p-1),8,:black)) ] )
    end
    Plots.annotate!(pltProb,[ ((propStats.T/20),N*1.07,("Z-mer",8,:black)) ] )
    return(pltProb)
end

function depthHistogramAnim(sims,par;colPalette = :cool,fps=60,secs=length(sims.Time)/fps,
    videoName="animation/3D/depthHistogramAnim.mp4",minimum = true,lLl = length(par.layerR[]),
    colors = palette(colPalette,lLl))
    
    N = size(sims.P,2)
    layerLims = [par.layerR[];par.obj.radius]
    radius=par.obj.radius

    I = round.(Int,range(0,1,round(Int,secs*fps)+1)[Not(1)] .* length(sims.Time))
    if minimum
        minDepthInfo = getTimeMinMatrix(sims,par.obj.radius,Vars=sims.minDepth)
        video = @animate for i in ProgressBar(I)
            depthHistogramPlot(minDepthInfo.sortedValues[i,:],i,layerLims,sims,par,radius=radius,
                N=N,lLl=lLl,colors=colors)
        end
    else
        posDistance = [L2(p) for p in sims.P]
        video = @animate for i in ProgressBar(I)
            depthHistogramPlot(posDistance[i,:],i,layerLims,sims,par,radius=radius,
                N=N,lLl=lLl,colors=colors,title="Radicle Distance from Center Histogram")
        end
    end
    mp4(video,videoName,fps=fps)
    return nothing
end

function depthHistogramPlot(values,i,layerLims,sims,par;radius=par.obj.radius,N = size(sims.P,2),lLl = length(par.layerR[]),
    colors = palette(:cool,lLl),title="Radicle Closest distance from Center Histogram")

    # standard histogram plots
    plt = Plots.histogram(values, label = "Histogram",bins=max(1,Int(floor(maximum(values) - minimum(values)))),
    xlabel = "Distance (nm)",ylabel = "Number of samples",title=title,xlims = (0,radius),ylims = (0,N))
    
    # Write time and zmer length
    Plots.annotate!(plt,[ (radius/2,N/1.04,(@sprintf("Time: %.7f sec",sims.Time[i]),10,:black)) ] )
    Plots.annotate!(plt,[ (radius/2,N/1.1,("Oligomer: "*string(sims.Zmer[i,1]),10,:black)) ] )
    
    # Plot layers with colors
    for l in 1:lLl
        BoxShape(plt,layerLims[[l,l+1]],[0,N];bw= 1,lc=plot_color(colors[l], 0.2),c = plot_color(colors[l], 0.2),label="Layer-"*string(l))
    end
    
    return(plt)
end

function BoxShape(plt,w,h;bw= 3,lc=:black,c = plot_color(:lightgreen, 0.1),label="")
    Plots.plot!(plt,Shape([(w[2],h[2]),(w[2],h[1]),(w[1],h[1]),(w[1],h[2])]),label=label,linewidth = bw,lc=lc,c=c)
    return(plt)
end

function AddVelocityBoxPlot(plt,Q,rmsv,Ermsv,w;linewidth = 1,rmsw= 3.0,title = "box",
    label = "label",ls=:dash,lc=:blue,lb=:black,fill=:white)
    # Add verticle line
    midw = mean(w)
    Plots.plot!(plt,[midw,midw],[Q[1];Q[2]] ,linewidth=linewidth,ls=ls,lc=lb,label="")
    Plots.plot!(plt,[midw,midw],[Q[3];Q[4]] ,linewidth=linewidth,ls=ls,lc=lb,label="")
    # Add horizontal Quantile Lines
    Plots.scatter!(plt,[midw],[Q[1]],m=:utriangle,markersize=rmsw,color=fill,label="")
    Plots.scatter!(plt,[midw],[Q[4]],m=:dtriangle,markersize=rmsw,color=fill,label="")
    # Plot empty box
    BoxShape(plt,w,[Q[2],Q[3]],bw=linewidth,lc=lb,c=fill,label =label)
    # Add horizontal Expected RMS line
    Plots.scatter!(plt,[midw],[Ermsv] ,markersize=rmsw,m=:cross,color=lc,label="")
    Plots.scatter!(plt,[midw],[rmsv] ,markersize=rmsw,m=:xcross,color=lb,label="")
    return(plt)
end

# This takes very long to plot due to the number of points. 
# return plts : a vector of plots.
function velocityScatter(sims)
    L2CMat = [L2Distance([0.0,0,0],vec(rad_p)) for rad_p in sims.P]
    VtMat  = [L2Distance([0.0,0,0],vec(v)) for v in sims.V]
    uZ = unique(sims.Zmer)
    plts = []
    for z in uZ
        TI = findall(sims.Zmer .== z)
        L2center = L2CMat[TI]
        VoverT = VtMat[TI]
        plt = Plots.scatter(VoverT,L2center,xlab="Velocity",ylab="Distance from center")
        for r in layerR[]
            Plots.plot!(plt,[minimum(VoverT),maximum(VoverT)],[r,r])
        end
        push!(plts,plt)
    end
    return plts
end

function velocityBoxPlot(sims,par,propStats,reactionTemp,Tg₀,τ;
    center=[0.0,0,0],bar_width = 0.8,qs = [0.05,0.25,0.75,0.95],Tmon=106,
    xlabel="Zmer length",ylabel = "velocity",title="Velocity Boxplot",leg_title ="Stats",lLl = length(par.layerR[]),
    colors = palette(:cool,lLl) )
    
    # get readius and speed from simulated plot and velocity
    L2CMat = [L2Distance(center,vec(rad_p)) for rad_p in sims.P]
    VtMat  = [L2Distance([0.0,0,0],vec(v)) for v in sims.V]

    # Initialize limits and lzmerl from samples.
    layerLims = [par.layerR[];par.obj.radius]
    uZ = collect(minimum(sims.Zmer)-1 .+ (1:length(propStats.wps)))

    # Calculating expected velocity # For Tg change with Wp : use (reactionTemp .- updateTg($Wp,Tg₀,Tmon))
    deltaDs = permutedims(hcat([[10^logD(wp,t,unit="nm") for t in (reactionTemp .- Tg₀)] for wp in propStats.wps]...))
    Dmat = deltaDs ./ (unique(uZ) .^ (0.5 .+ 1.75 .* propStats.wps ))
    Ermsv = [sqrt(6 * eD / τ[]) for eD in Dmat]

    width_i = collect(1:(lLl+1)) ./ (lLl+2)
    width_delta = 1/(lLl+2)

    bp = Plots.plot(
        xlabel=xlabel,ylabel = ylabel,title=title,leg_title =leg_title, dpi=300, yaxis=:log10,
        xlims =(minimum(uZ)-1,maximum(uZ)+1),ylims = (quantile(reduce(vcat,VtMat),qs[1]/length(par.layerR[])),maximum(VtMat)),
        legend=:outerright,xticks = uZ
    )
    # Group by zmer length
    for z in unique(uZ)
        TI = findall(sims.Zmer .== z)
        L2center = L2CMat[TI]
        VoverT = VtMat[TI]
        # group by layers
        for layer in 1:lLl
            l_lower = layerLims[layer]
            l_upper = layerLims[layer+1]
            LI = findall(l_lower.< L2center .<= l_upper)
            if !isempty(LI)
                V2CL = VoverT[LI]
                rmsv = sqrt(mean(V2CL.^2))
                w1 = z - 0.5 + width_i[layer]+ width_delta*(1-bar_width)/2
                w2 = z - 0.5 + width_i[layer+1]- width_delta*(1-bar_width)/2
                Q = quantile(V2CL,qs)
                AddVelocityBoxPlot(
                    bp,Q,rmsv,Ermsv[z-minimum(uZ)+1,layer],[w1,w2];
                    linewidth = 1,rmsw= 2,label = "",
                    ls=:solid,lc=:black,lb=:black,fill=colors[layer]
                )
            end
        end
    end

    for l in 1:lLl
        BoxShape(bp,(minimum(uZ) .- [100,99]),(maximum(VtMat)*10 .+ [1,2]);
        bw= 1,c = colors[l],label="Layer-"*string(l)*" IQR")
    end
    Plots.scatter!(bp,[],[],label="Expected RMS",m=:cross,color=:black)
    Plots.scatter!(bp,[],[],label="Simulate RMS",m=:xcross,color=:black)
    Plots.scatter!(bp,[],[],label=string(round(Int,last(qs)*100))*"%-quantile",m=:dtriangle,color=:white,linewidth =1)
    Plots.scatter!(bp,[],[],label=string(round(Int,first(qs)*100))*"%-quantile",m=:utriangle,color=:white,linewidth =1)
    return(bp)
end
# function distanceAnim(sims,par,j,propTime,lPropl,lZmerl;videoSec = 20,fps = 60)
#     I = [1;round.(Int,range(0,1,videoSec*fps)[Not(1)] .* length(sims.Time))]
#     video = @animate for i in ProgressBar(I)
#         distancePlot(sims,par,j,propTime,lPropl,lZmerl,i=i)
#     end
#     return(video)
# end

# function distancePlot(
#     sims,par,j,propTime,lPropl,lZmerl;i=length(sims.Time),center = [0.0,0,0],lLl = length(par.layerR[]),
#     colors = palette(:cool, lLl)) # ,colors = range(HSV(-120,1,1), stop=HSV(-360,1,1), length=lLl)

#     layerLims = [par.layerR[];par.obj.radius]
#     L2vj=[L2Distance(center,vec(rad_p)) for rad_p in sims.P[:,j]]
#     Height = [0,par.obj.radius]
#     plt = Plots.plot(
#         xlim=(0,last(sims.Time)),ylim=(0,par.obj.radius*1.01),
#         legend=:outerright , xlabel = "time (second)",ylabel = "Distance from center (nm)",
#         title = "Radical distance from center",size = (700,600), xticks=range(0,round.(last(sims.Time),sigdigits=2),6)[Not(1)]
#     )
#     # Add layer colors
#     for l in 1:lLl
#         BoxShape(plt,[0,last(sims.Time)],layerLims[[l,l+1]];bw= 1,lc=plot_color(colors[l], 0.2),c = plot_color(colors[l], 0.2),label="Layer-"*string(l))
#     end
#     # Add Zmer Length
#     for p in 0:(lPropl-1)
#         t_i = propTime * (p+1)
#         Plots.plot!(plt,[t_i,t_i],Height,lc=plot_color(:grey,0.5),label="")
#         Plots.annotate!(plt,[ (t_i-0.5*propTime,Height[2]*1.01,(string(lZmerl+p),8,:black)) ] )
#     end
#     Plots.annotate!(plt,[ (propTime*(lPropl*1.07),Height[2]*1.01,("Z-mer",8,:black)) ] )

#     Plots.plot!(plt,sims.Time[1:i],L2vj[1:i],lc=:black,label="radicle")
#     return(plt)
# end

function addToDict(d, wpInit,totalLayer, glossy,value)
    if !(string(wpInit) in keys(d))
        d[string(wpInit)] = Dict()
    end
    if !(string(totalLayer) in keys(d[string(wpInit)]))
        d[string(wpInit)][string(totalLayer)] = Dict()
    end
    d[string(wpInit)][string(totalLayer)][string(glossy)] = value
end

function check_path(path::AbstractString)
    directory = join(split(path, "/")[begin:end-1],"/")
    if (directory != "") && (!isdir(directory))  # Check if the directory does not exist
        mkpath(directory) 
    end
    return path
end

function init_jld(filename)
    if isfile(filename)
        return load(filename)
    else 
        return Dict()
    end
end

function save_jld(filename,object)
    check_path(filename)
    save(filename,object)
end

function load_jld(filename)
    return load(filename)
end

function wdir(folder)
    if !isempty(folder) && !ispath(folder)
        println("Create directory ("*folder*")!")
        mkpath(folder)
    end
    return(folder)
end

function writeSimsData(sims;folder="",overwrite=false)
    cur = ""
    for fold in split(folder, "/")
        cur *= (fold * "/")
        wdir(cur)
    end
    cur = ( length(cur)>1 ? cur : "" )
    if overwrite || !ispath(cur*"data1.csv")
        N = size(sims.P,2)
        for i in 1:N
            Ps = reduce(hcat,[vec(p) for p in sims.P[:,i]])
            Vs = reduce(hcat,[vec(v) for v in sims.V[:,i]])
            df = DataFrame(
                x=Ps[1,:],y=Ps[2,:],z=Ps[3,:],vx=Vs[1,:],vy=Vs[2,:],vz=Vs[3,:],
                Zmer=sims.Zmer[:,i],minL=sims.minL[:,i],Time=sims.Time)
            CSV.write(cur*"data"*string(i)*".csv", df)
        end
    else
        println("The files ("*cur*"data1.csv) exist. Set overwrite=true, if want it to be overwrite.")
    end
end


