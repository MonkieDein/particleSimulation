include("chem.jl")
using Plots
using Colors
using CSV

function sampleDepthPlot(sims,lLl,lPropl,propTime,lZmerl)
    lTl = length(sims.Time)
    N = size(sims.P,2)
    Nreach = zeros(Int,size(sims.Time))
    Layers = ones(Int,(lTl,N)) .+ lLl
    sortLayers = ones(Int,(lTl,N)) .+ lLl
    for i in 1:lTl
        Layers[i,:] = min.(Layers[max(i-1,1),:],sims.minL[i,:])
        sortLayers[i,:] = sort(copy(Layers[i,:]))
        Nreach[i] = sum(Layers[i,:] .== 1)
    end
    # Plot Heat Map
    pltProb = heatmap(sims.Time, 1:N,transpose(sortLayers),color = palette(:cool, lLl),
    title="Distribution of Deepest Layer reached",legend=:outerright,
    xlabel="Time", ylabel="No of samples",ylims=(0,N*1.1),colorbar=false,leg_title ="Layer" , dpi=300
    )
    # Plot Label
    for (i,c) in enumerate(palette(:cool, lLl))
        scatter!(pltProb,[],[],color=c,label=string(i))
    end
    # Add Zmer Lenght
    for p in 0:(lPropl-1)
        t_i = propTime * (p+1)
        plot!(pltProb,[t_i,t_i],[0,N],lc=:black,label="")
        annotate!(pltProb,[ (t_i-0.5*propTime,105,(string(lZmerl+p),8,:black)) ] )
    end
    annotate!(pltProb,[ (propTime*(lPropl+0.5),105,("Z-mer",8,:black)) ] )
    return(pltProb)
end

function BoxShape(plt,w,h;bw= 3,lc=:black,c = plot_color(:lightgreen, 0.1),label="")
    plot!(plt,Shape([(w[2],h[2]),(w[2],h[1]),(w[1],h[1]),(w[1],h[2])]),label=label,linewidth = bw,lc=lc,c=c)
    return(plt)
end

function AddVelocityBoxPlot(plt,Q,rmsv,Ermsv,w;linewidth = 1,rmsw= 3.0,title = "box",
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
        xlabel=xlabel,ylabel = ylabel,title=title,leg_title =leg_title, dpi=300,
        xlims =(minimum(uZ)-1,maximum(uZ)+1),legend=:outerright,xticks = uZ
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
            V2CL = VoverT[LI]
            rmsv = sqrt(mean(V2CL.^2))
            w1 = z - 0.5 + width_i[layer]+ width_delta*(1-bar_width)/2
            w2 = z - 0.5 + width_i[layer+1]- width_delta*(1-bar_width)/2
            Q = quantile(V2CL,qs)
            AddVelocityBoxPlot(
                bp,Q,rmsv,Ermsv[z-minimum(uZ)+1,layer],[w1,w2];
                linewidth = 1,rmsw= 2,label = ( z==first(uZ) ? "Layer-"*string(layer)*" IQR" : ""),
                ls=:solid,lc=:red,lb=:black,fill=colors[layer]
            )
        end
    end
    scatter!(bp,[],[],label="Expected RMS",m=:cross,color=:red)
    scatter!(bp,[],[],label="Simulate RMS",m=:xcross,color=:black)
    scatter!(bp,[],[],label=string(round(Int,last(qs)*100))*"%-quantile",m=:dtriangle,color=:white,linewidth =1)
    scatter!(bp,[],[],label=string(round(Int,first(qs)*100))*"%-quantile",m=:utriangle,color=:white,linewidth =1)
    return(bp)
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


