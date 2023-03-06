include("chem.jl")
using Plots
using Colors
using CSV


function simpleDepthPlot(
    sims,lLl;lTl= length(sims.Time),N= size(sims.P,2),
    colors = palette(:cool, lLl),title="", xlabel="", ylabel="")
    Layers = ones(Int,(lTl,N)) .+ lLl
    sortLayers = ones(Int,(lTl,N)) .+ lLl
    for i in 1:lTl
        Layers[i,:] = min.(Layers[max(i-1,1),:],sims.minL[i,:])
        sortLayers[i,:] = sort(copy(Layers[i,:]))
    end
    existLayers = collect(minimum(sortLayers):maximum(sortLayers))
    # Plot Heat Map
    pltProb = heatmap(
        sims.Time, 1:N,transpose(sortLayers),color = colors[existLayers] ,
        colorbar=false, dpi=300,title=title,xlabel=xlabel,ylabel=ylabel)
    return(pltProb)
end

function sampleDepthPlot(sims,lLl,lPropl,propTime,lZmerl)
    lTl = length(sims.Time)
    N = size(sims.P,2)
    colors = palette(:cool, lLl)

    # Plot Heat Map
    pltProb = simpleDepthPlot(sims,lLl,lTl=lTl,N=N,colors=colors,
    title="Distribution of Deepest Layer reached",
    xlabel="Time", ylabel="No of samples")
    
    plot!(pltProb,ylims=(0,N*1.1),legend=:outerright,leg_title ="Layer")
    # Plot Label
    for (i,c) in enumerate(colors)
        scatter!(pltProb,[],[],color=c,label=string(i))
    end
    # Add Zmer Length
    for p in 0:(lPropl-1)
        t_i = propTime * (p+1)
        plot!(pltProb,[t_i,t_i],[0,N],lc=:black,label="")
        annotate!(pltProb,[ (t_i-0.5*propTime,N*1.05,(string(lZmerl+p),8,:black)) ] )
    end
    annotate!(pltProb,[ (propTime*(lPropl+0.5),N*1.05,("Z-mer",8,:black)) ] )
    return(pltProb)
end


function BoxShape(plt,w,h;bw= 3,lc=:black,c = plot_color(:lightgreen, 0.1),label="")
    plot!(plt,Shape([(w[2],h[2]),(w[2],h[1]),(w[1],h[1]),(w[1],h[2])]),label=label,linewidth = bw,lc=lc,c=c)
    return(plt)
end

function AddVelocityBoxPlot(plt,Q,rmsv,Ermsv,w;linewidth = 1,rmsw= 3.0,
    label = "label",ls=:solid,lc=:blue,lb=:black,fill=:white)
    # Add verticle line
    midw = mean(w)
    plot!(plt,[midw,midw],[Q[1];Q[2]] ,linewidth=linewidth,ls=:dash,lc=lb,label="")
    plot!(plt,[midw,midw],[Q[3];Q[4]] ,linewidth=linewidth,ls=:dash,lc=lb,label="")
    # Add horizontal Quantile Lines
    scatter!(plt,[midw],[Q[1]],m=:utriangle,markersize=rmsw,color=fill,label="")
    scatter!(plt,[midw],[Q[4]],m=:dtriangle,markersize=rmsw,color=fill,label="")
    # Plot empty box
    BoxShape(plt,w,[Q[2],Q[3]],bw=linewidth,lc=lb,c=fill,label =label)
    # Add horizontal Expected RMS line
    scatter!(plt,[midw],[Ermsv] ,markersize=rmsw,m=:cross,color=lc,label="")
    scatter!(plt,[midw],[rmsv] ,markersize=rmsw,m=:xcross,color=lb,label="")
    return(plt)
end

function velocityBoxPlot(sims,par;center=[0.0,0,0],bar_width = 0.8,qs = [0.05,0.25,0.75,0.95],
    xlabel="Zmer length",ylabel = "velocity",title="Velocity Boxplot",leg_title ="Stats")
    lLl = length(par.L.R)
    L2CMat = [L2Distance(center,vec(rad_p)) for rad_p in sims.P]
    VtMat  = [L2Distance([0.0,0,0],vec(v)) for v in sims.V]

    layerLims = [par.L.R;par.obj.radius]
    uZ = unique(sims.Zmer)
    # Dmat[zmer-zstart+1,layer] # empiricle expected RMS velocity
    Dmat = par.L.D' ./ (unique(uZ) .^ (0.5+1.75*par.Wp))
    Ermsv = [sqrt(6 * eD / τ) for eD in Dmat]

    colors = palette(:cool, lLl)
    width_i = collect(1:(lLl+1)) ./ (lLl+2)
    width_delta = 1/(lLl+2)

    bp = plot(
        xlabel=xlabel,ylabel = ylabel,title=title,leg_title =leg_title, yaxis=:log10,dpi=300,
        xlims =(minimum(uZ)-1,maximum(uZ)+1),xticks = uZ,legend=:outerright
    )#,ylims = (1,maximum(VtMat)),yticks = 10 .^ (1:ceil(log(maximum(VtMat))/log(10)))
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
                    ls=:solid,lc=:red,lb=:black,fill=colors[layer]
                )
            end
        end
    end
    # Use Number beyond xlims and ylims
    for layer in 1:lLl
        BoxShape(bp,[maximum(uZ)+99;maximum(uZ)+100],[maximum(VtMat);maximum(VtMat)+1];bw= 1,lc=:black,c = colors[layer],label="Layer-"*string(layer)*" IQR")
    end
    scatter!(bp,[],[],label="Expected RMS",m=:cross,color=:red)
    scatter!(bp,[],[],label="Simulate RMS",m=:xcross,color=:black)
    scatter!(bp,[],[],label=string(round(Int,last(qs)*100))*"%-quantile",m=:dtriangle,color=:white,linewidth =1)
    scatter!(bp,[],[],label=string(round(Int,first(qs)*100))*"%-quantile",m=:utriangle,color=:white,linewidth =1)
    return(bp)
end
function distanceAnim(sims,par,j,propTime,lPropl;videoSec = 20,fps = 60)
    I = [1;round.(Int,range(0,1,videoSec*fps)[Not(1)] .* length(sims.Time))]
    distancePlot(sims,par,j,propTime,lPropl)

    video = @animate for i in ProgressBar(I)
        distancePlot(sims,par,j,propTime,lPropl,i=i)
    end
    return(video)
end

function distancePlot(
    sims,par,j,propTime,lPropl;i=length(sims.Time),center = [0.0,0,0])

    lLl = length(par.L.R)
    colors = palette(:cool, lLl)
    layerLims = [par.L.R;par.obj.radius]
    L2vj=[L2Distance(center,vec(rad_p)) for rad_p in sims.P[:,j]]
    Height = [0,par.obj.radius]
    plt = plot(
        xlim=(0,last(sims.Time)),ylim=(0,par.obj.radius*1.01),
        legend=:outerright , xlabel = "time (second)",ylabel = "Distance from center (nm)",
        title = "Radical distance from center",size = (700,600), xticks=range(0,round.(last(sims.Time),sigdigits=2),6)[Not(1)]
    )
    # Add layer colors
    for l in 1:lLl
        BoxShape(plt,[0,last(sims.Time)],layerLims[[l,l+1]];bw= 1,lc=plot_color(colors[l], 0.2),c = plot_color(colors[l], 0.2),label="Layer-"*string(l))
    end
    # Add Zmer Length
    for p in 0:(lPropl-1)
        t_i = propTime * (p+1)
        plot!(plt,[t_i,t_i],Height,lc=plot_color(:grey,0.5),label="")
        annotate!(plt,[ (t_i-0.5*propTime,Height[2]*1.01,(string(lZmerl+p),8,:black)) ] )
    end
    annotate!(plt,[ (propTime*(lPropl+0.5),Height[2]*1.01,("Z-mer",8,:black)) ] )

    plot!(plt,sims.Time[1:i],L2vj[1:i],lc=:black,label="radicle")
    return(plt)
end

function wdir(folder)
    if !isempty(folder) && !ispath(folder)
        println("Create directory ("*folder*")!")
        mkdir(folder)
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


