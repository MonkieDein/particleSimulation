include("3Dbasics.jl")

using DataFrames
using GLMakie
using Colors

# Region is treated as a constant
Region = DataFrame(
    C1 = [ -4.428 ; 26.0 ; 159.0 ; -13.7 ],
    C2 = [ 1.842 ; 37.0 ; 170.0 ; 0.500 ],
    C3 = [ 0 ; 0.0797 ; 0.3664 ; 0 ],
    C4 = [ 8.12e-3 ; 0 ; 0 ; 0 ]
)

# logD = ( A - B ⋅ wₚ ) where A = C₁ + ΔT C₃ , and B = C₂ - ΔT C₄
# ΔT = Tᵣₓₙ - Tg : Tg = 70
function logD( wₚ , ΔT; unit = "cm", Region = Region)   
    scale = Dict("m" => -4,"cm" => 0,"nm" => 14)
    A = Region.C1 .+ ΔT .* Region.C3  
    B = Region.C2 .- ΔT .* Region.C4 
    # Intersection of wₚ : (aᵢ - aᵢ₊₁)/(bᵢ - bᵢ₊₁) or (aᵢ₊₁ - aᵢ)/(bᵢ₊₁ - bᵢ)
    RegionMin = [0 ; (A[[1,2,3]] .- A[[2,3,4]]) ./ (B[[1,2,3]] .- B[[2,3,4]])]
    RegionIndex = searchsortedlast(RegionMin,wₚ)
    return( ( A[RegionIndex] - B[RegionIndex] * wₚ ) + scale[unit])
end

mutable struct mLparticle
    obj :: Round # particle itself as a Round object
    Wp :: Observable{Float64} # Weightage of polymer
    layerR :: Observable{Vector{Float64}} # Layers Radius start
    layerΔT :: Observable{Vector{Float64}} # Layers Delta T
    layerD :: Observable{Vector{Float64}} # Layers Diffusion coefficient
    function mLparticle(R,Wp,layerR,layerΔT,layerD;x=0.0,y=0.0,z=0.0)
        new(Round(x,y,z,R),Wp,layerR,layerΔT,layerD)
    end
end

function autoUpdateRadVelocity(D::Observable{Float64},obj::Round,σx::Observable{Float64},τ::Observable{Float64})
    on(D) do Dval
        newσx = sqrt( 2 * Dval / τ[]);
        obj.v = coord( vec(obj.v) .* (newσx / σx[]) ); # scale v multiplying (new_σx/old_σx) 
        σx[] = newσx;
    end

    on(τ) do Tval
        newσx = sqrt( 2 * D[] / Tval);
        obj.v = coord( vec(obj.v) .* (newσx / σx[]) ); # scale v multiplying (new_σx/old_σx) 
        σx[] = newσx;
    end
end

mutable struct Radicle
    obj :: Round # radicle itself as a Round object
    D :: Observable{Float64} # diffusion coefficient
    zmer :: Observable{Int} # zipper mer length
    τ :: Observable{Float64} # Global Time step
    σx :: Observable{Float64} # Mean square displacement of each axis
    l :: Observable{Int} # layer index
    function Radicle(x,y,z,par;r=0,zmer = Observable(1),τ = Observable(0.01))
        l = Observable( searchsortedlast(par.layerR[], L2Distance( [x,y,z],vec(par.obj.p) ) ) )
        D = @lift( $(par.layerD)[$l]/($zmer^(0.5+1.75*$(par.Wp) )) ) 
        # lift is not used in "σx" because "obj.v" is depend on the previous "σx". 
        σx = Observable(sqrt( 2 * D[] / τ[]))
        obj = Round([x,y,z],r,V=[σx[],0.0,0.0])
        autoUpdateRadVelocity(D,obj,σx,τ)
        new(obj,D,zmer,τ,σx,l)
    end
    function Radicle(P,par;r=0,l = Observable(1),zmer = Observable(1),τ = Observable(0.01))
        l = Observable( searchsortedlast(par.layerR[], L2Distance( P,vec(par.obj.p) ) ) )
        D = @lift( $(par.layerD)[$l]/($zmer^(0.5+1.75*$(par.Wp) )) ) 
        # lift is not used in "σx" because "obj.v" is depend on the previous "σx". 
        # need an ObservableFunction to auto update radicle "velocity" and "σx".
        σx = Observable(sqrt( 2 * D[] / τ[]))
        obj = Round(P,r,V=[σx[],0.0,0.0])
        autoUpdateRadVelocity(D,obj,σx,τ)
        new(obj,D,zmer,τ,σx,l)
    end
end

function °C2K(T)
    return(T .+ 273.15)
end

function K2°C(T)
    return(T .- 273.15)
end

function updateTg(wp,tgseed,tmon)
    return K2°C( 1 ./ ( (wp ./ °C2K(tgseed)) .+ ((1-wp) ./ °C2K(tmon)) ))
end


function LayerTransitionUpdate(particle::mLparticle,radicle::Radicle,l::Int)
    d = normalize( vec(particle.obj.p) .- vec(radicle.obj.p) )
    v = normalize( vec(radicle.obj.v) ) 
    radicle.l[] = max( 1,l - (dot(d,v) > 0) ) 
end

# If postTransition is True ⟹ min is zero and we already update the transition
function multiLbounceStepUpdate(particle::mLparticle,radicle::Radicle,Δt;postTransition = false,wallCollideDo = "reflect")  
    tempΔt = Δt
    Center = vec(particle.obj.p)
    InitPoint = vec(radicle.obj.p)
    minDepth = L2Distance(Center,InitPoint)
    while tempΔt > 0
        # I = Layer's intersection within time portion [0,1). 1 ⟹ no intersection.
        I = reduce(vcat,[intersections(vec(particle.obj.p), R, radicle.obj, tempΔt,collide=false) for R in particle.layerR[]])
        # l : the index of (1) the smallest element or (2) the second smallest element
        l = (1 + partialsortperm( I, (1 + postTransition) )) ÷ 2
        transition_time = I[partialsortperm( I, (1 + postTransition) )]
        if transition_time == 0        # Transitional update
            # println("Transition update",tempΔt)
            LayerTransitionUpdate(particle,radicle,l)
            Depth = multiLbounceStepUpdate(particle,radicle,tempΔt,postTransition = true,wallCollideDo = wallCollideDo)
            minDepth = min(Depth,minDepth)
            tempΔt = 0
        elseif transition_time < 1     # Intersect with another Layer 
            # println("Intersect with another Layer ",tempΔt ," position ",vec(radicle.obj.p)," trans time ", transition_time)
            preUpdatePosition = vec(radicle.obj.p)
            preL2 = L2(radicle.obj.p)
            updateMotion(radicle.obj,transition_time* tempΔt)
            tempΔt -= transition_time * tempΔt
            # if updateMotion is way too small, no positional change
            # handle 1e16 loop -> make a tiny jump instead
            if ( ((sum(abs.(vec(radicle.obj.p) .- preUpdatePosition)) == 0) || (L2(radicle.obj.p) == preL2 ) ) && 
                (transition_time < 1e-10) )
                # println("eps intersect",tempΔt)
                jumptime = min(1/(L2(radicle.obj.v)*1e8),tempΔt/1e10)
                updateMotion(radicle.obj,jumptime)
                tempΔt -= jumptime
            end
            minDepth = min(minDepth,closestDistance(Center,InitPoint,vec(radicle.obj.p)))
            
            # FP addition is not associative, transition update is necessary to avoid surpass layer
            LayerTransitionUpdate(particle,radicle,l)
        else                           # Collide on particle wall
            colission = collisionNreflection(particle.obj,radicle.obj,tempΔt)
            updateMotion(radicle.obj,colission.time * tempΔt)
            minDepth = min(minDepth,closestDistance(Center,InitPoint,vec(radicle.obj.p)))
            if (colission.time < 1)
                if wallCollideDo == "reflect"
                    radicle.obj.v = colission.reflection
                elseif wallCollideDo == "rand_dir"
                    while L2Distance(vec(particle.obj.p),vec(radicle.obj.p) .+ (vec(radicle.obj.v) ./ (L2(radicle.obj.v) * 1e8)) ) ≥ (particle.obj.radius - radicle.obj.radius)
                        # println("rand_dir on wall collision",tempΔt)
                        radicle.obj.v = random_direction(radicle.obj,size=radicle.σx[])
                    end
                elseif wallCollideDo == "random"
                    while L2Distance(vec(particle.obj.p),vec(radicle.obj.p) .+ (vec(radicle.obj.v) ./ (L2(radicle.obj.v) * 1e8)) ) ≥ (particle.obj.radius - radicle.obj.radius)
                        # println("random on wall collision",tempΔt)
                        radicle.obj.v = random_velocity(radicle.obj,radicle.σx[])
                    end
                end
            end
            tempΔt -= colission.time * tempΔt
        end
        postTransition = false # reset postTransition
        Center = vec(particle.obj.p)
        InitPoint = vec(radicle.obj.p)
    end
    return(minDepth)
end

# Propagation Time for Monomer
function propagationTimeInterval(Wp,TempInC;Mw = 142.2  )
    R = 8.314;                                              # Gas Constant (J/(K * mol))
    Mw = Mw                                                # Molecule weight (g/mol) ! 
    M = ((1-Wp)/Mw)*1000;                                   # Monomer concentration (mol/g) !
    k_p = (2.673e6)*exp(-22.36e3/(R * °C2K(TempInC)));      # Propagation Coefficient ( L/(mol*s) )
    propTime = 1/(k_p * M);                                 # Propagation Time
    return(propTime)
end

function wpPiecewiseLinear(wpInit , endTime, curTime ; wpEnd = 1.0, startTime=0.0)
    if curTime < startTime
        return wpInit
    elseif curTime > endTime
        return wpEnd
    else
        return wpInit + (wpEnd - wpInit) * (curTime - startTime) / (endTime - startTime)
    end
end

function getPropStats(wpInit,reactionTemp,endTime,nProp; wpEnd = 1.0, startTime=0.0,currentTime = 0.0)
    times = []
    wps = []
    wp = wpPiecewiseLinear(wpInit,endTime,currentTime,wpEnd=wpEnd,startTime=startTime)
    for prop in 1:nProp
        push!(wps,wp)
        propTime = propagationTimeInterval(wp,reactionTemp)
        push!(times,propTime)
        currentTime = currentTime + propTime
        wp = wpPiecewiseLinear(wpInit,endTime,currentTime,wpEnd=wpEnd,startTime=startTime)
    end
    return (wps=wps,times=times,T=currentTime)
end

#TODO: need to figure out how to get total time elapsed
function simulate(propTime::Observable{Float64},Wp::Observable{Float64},
    zmer::Observable{Int},Time,par::mLparticle,Rad::Vector{Radicle},
    wpInit::Float64,linReactEndTime::Float64;
    wpEnd=1.0,startTime=0.0,
    wallCollideDo = "reflect", eachStepDo ="random")

    # Stats of Simulation 
    N = length(Rad)
    τs = diff(Time)
    lTl = length(Time)                        
    pos = Array{coord}(undef, (lTl,N))      
    vel = Array{coord}(undef, (lTl,N))     
    zmers = zeros(Int,(lTl,N))
    minDepth = zeros(Float64,(lTl,N)) .+ par.obj.radius
    tnextP = propTime[]

    for (i,t) in ProgressBar(enumerate(Time))
        # ------ Collect Statistics -----
        for (j,rad) in enumerate(Rad)
            if eachStepDo == "random"
                random_velocity(rad.obj,rad.σx[])
            elseif eachStepDo == "rand_dir"
                random_direction(rad.obj,size=rad.σx[])
            end
            pos[i,j] = rad.obj.p 
            vel[i,j] = rad.obj.v 
            zmers[i,j] = rad.zmer[] 
            minDepth[i,j] = min( minDepth[i,j], L2Distance(vec(par.obj.p),vec(rad.obj.p)) )
        end    
        # ------ Reach End Time, Finish. -----
        if i == lTl
            break
        end
        # ------ Taking Time Step Update -----
        tempτ = τs[i]
        while (tnextP ≤ t  + tempτ) 
            prePropTime = tnextP - t # remaining time before propagation
            # # println("time step: ",tempτ,"\t ",prePropTime)
            # Step Update remaining time before propagation
            for (j,rad) in enumerate(Rad)
                # ----- udpate position before propagation -----
                curL = multiLbounceStepUpdate(par,rad,prePropTime,wallCollideDo=wallCollideDo)
                minDepth[i+1,j] = min( minDepth[i+1,j] ,curL)
            end
            tempτ -= prePropTime
            # ----- radicle undergo propagation -----
            Wp[] = wpPiecewiseLinear(wpInit,linReactEndTime,tnextP,wpEnd=wpEnd,startTime=startTime)
            zmer[] += 1
            tnextP += propTime[]
        end
        # ----- udpate position after propagation -----
        for (j,rad) in enumerate(Rad)
            curL = multiLbounceStepUpdate(par,rad,tempτ,wallCollideDo=wallCollideDo)
            minDepth[i+1,j] = min( minDepth[i+1,j] ,curL)
        end
    end
    
    return(P=pos,V=vel,Zmer=zmers,Time=Time,minDepth=minDepth,
    minL=[searchsortedlast(par.layerR[], d ) for d in minDepth])
end


function MaxStepSize(D,Wp,zmer,τ;sd=1,confident=(1-cdf(Distributions.Normal(),-sd)*2))
    Dmax = maximum(D)/(zmer^(0.5+1.75*Wp));
    L2_std = sqrt( 6 * Dmax * τ ); 
    q = (1 .- confident) ./ 2
    Z = abs.(Distributions.quantile.(Distributions.Normal(0,L2_std), q))
    return(Z)
end

function MinTimeForStepsize(stepsize,D,Wp,zmer;
    sd=1,confident=(1-cdf(Distributions.Normal(),-sd)*2))
    τ = 1
    Dmax = maximum(D)/(zmer^(0.5+1.75*Wp));
    L2_std = sqrt( 6 * Dmax * τ ); 
    q = (1 .- confident) ./ 2
    Z = abs.(Distributions.quantile.(Distributions.Normal(0,L2_std), q))
    return( (( stepsize ./ Z' ) .^ 2) )
end
