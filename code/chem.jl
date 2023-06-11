include("3Dbasics.jl")
using DataFrames

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
    Wp :: Float64 # Weightage of polymer
    L :: DataFrame # Layers rstart and Diffusion coefficient
    function mLparticle(R,Wp,L;x=0.0,y=0.0,z=0.0,D = [10^logD(Wp,t,unit="nm") for t in L.ΔT])
        L.D = D
        new(Round(x,y,z,R),Wp,L)
    end
end

mutable struct Radicle
    obj :: Round # radicle itself as a Round object
    D :: Float64 # diffusion coefficient
    zmer :: Int # zipper mer length
    τ :: Float64 # Global Time step
    σx :: Float64 # Mean square displacement of each axis
    l :: Int # layer index
    function Radicle(x,y,z;D = 1e7,r=0,vx=1,vy=0,vz=0,σx = 1.0,l = 1,zmer = 1,τ = 0.01)
        new(Round([x,y,z],r,V=[vx,vy,vz]),D,zmer,τ,σx,l)
    end
    function Radicle(P;D = 1e7, r=0,V=[1.0,0,0],σx = 1.0,l = 1,zmer = 1,τ = 0.01)
        new(Round(P,r,V=V),D,zmer,τ,σx,l)
    end
end

function °C2K(T)
    return(T + 273.15)
end

function K2°C(T)
    return(T - 273.15)
end


function updateRadicle(particle::mLparticle,radicle::Radicle )
    D = particle.L.D[radicle.l];
    radicle.D = D/(radicle.zmer^(0.5+1.75*particle.Wp));
    newσx = sqrt( 2 * radicle.D / radicle.τ);
    radicle.obj.v = coord( vec(radicle.obj.v) .* (newσx / radicle.σx) ); # scale v multiplying (new_σx/old_σx) 
    radicle.σx = newσx;
end

function LayerTransitionUpdate(particle::mLparticle,radicle::Radicle,l::Int)
    d = normalize( vec(particle.obj.p) .- vec(radicle.obj.p) )
    v = normalize( vec(radicle.obj.v) ) 
    radicle.l = max( 1,l - (dot(d,v) > 0) ) # If dot product > 0 means going toward the center
    updateRadicle(particle,radicle) # Update radicle new layer information
end

# If postTransition is True ⟹ min is zero and we already update the transition
function multiLbounceStepUpdate(particle::mLparticle,radicle::Radicle,Δt;postTransition = false)  
    tempΔt = Δt
    minL = radicle.l
    while tempΔt > 0
        # I = Layer's intersection within time portion [0,1). 1 ⟹ no intersection.
        I = reduce(vcat,[intersections(vec(particle.obj.p), R, radicle.obj, tempΔt,collide=false) for R in particle.L.R])
        # l : the index of (1) the smallest element or (2) the second smallest element
        l = (1 + partialsortperm( I, (1 + postTransition) )) ÷ 2
        transition_time = I[partialsortperm( I, (1 + postTransition) )]
        if transition_time == 0        # Transitional update
            LayerTransitionUpdate(particle,radicle,l)
            L = multiLbounceStepUpdate(particle,radicle,tempΔt,postTransition = true)
            minL = min(L,minL)
            tempΔt = 0
        elseif transition_time < 1     # Intersect with another Layer 
            preUpdatePosition = vec(radicle.obj.p)
            updateMotion(radicle.obj,transition_time* tempΔt)
            # if updateMotion is way too small, no positional change
            # handle 1e16 loop -> incremental position of highest velocity axis
            if sum(abs.(vec(radicle.obj.p) .- preUpdatePosition)) == 0 
                V = vec(radicle.obj.v)
                (v,i) = findmax(abs.(V))
                if V[i] > 0
                    preUpdatePosition[i] = nextfloat(preUpdatePosition[i])
                else
                    preUpdatePosition[i] = prevfloat(preUpdatePosition[i])
                end
                radicle.obj.p = coord(preUpdatePosition)
            end
            tempΔt -= transition_time * tempΔt
            # FP addition is not associative, transition update is necessary to avoid surpass layer
            LayerTransitionUpdate(particle,radicle,l)
        else                           # Collide on particle wall
            colission = collisionNreflection(particle.obj,radicle.obj,tempΔt)
            updateMotion(radicle.obj,colission.time * tempΔt)
            if (colission.time < 1)
                radicle.obj.v = colission.reflection
            end
            tempΔt -= colission.time * tempΔt
        end
        postTransition = false # reset postTransition
    end
    return(minL)
end


# Propagation Time for Monomer
function propagationTimeInterval(Wp,TempInC)
    R = 8.314;                                              # Gas Constant (J/(K * mol))
    Mw = 100                                                # Molecule weight (g/mol) ! 
    M = ((1-Wp)/Mw)*1000;                                   # Monomer concentration (mol/g) !
    k_p = (2.673e6)*exp(-22.36e3/(R * °C2K(TempInC)));      # Propagation Coefficient ( L/(mol*s) )
    propTime = 1/(k_p * M);                                 # Propagation Time
    return(propTime)
end

function simulate(τ::Float64,propTime::Float64,T::Float64,par::mLparticle,Rad::Vector{Radicle})
    # Stats of Simulation 
    Time = collect(0:τ:T)
    lTl = length(Time)                        
    N = length(Rad)
    pos = Array{coord}(undef, (lTl,N))      
    vel = Array{coord}(undef, (lTl,N))     
    zmers = zeros(Int,(lTl,N))
    minL = ones(Int,(lTl,N)) .+ length(par.L.D)
    tnextP = zeros(N) .+ propTime
    for (i,t) in ProgressBar(enumerate(Time))
        for (j,rad) in enumerate(Rad)
            random_velocity(rad.obj,rad.σx)
            # ------ Collect Statistics -----
            pos[i,j] = rad.obj.p 
            vel[i,j] = rad.obj.v 
            zmers[i,j] = rad.zmer 
            minL[i,j] = rad.l
            # ------  ----- ----- -----  -----
            # ------ Taking Time step update -----
            tempτ = τ
            while (tnextP[j] - t <= tempτ) 
                # ----- udpate position before propagation -----
                prePropTime = tnextP[j] - t # remaining time before propagation
                curL = multiLbounceStepUpdate(par,rad,prePropTime)
                minL[i,j] = min( minL[i,j] ,curL)
                tempτ -= prePropTime
                # ----- radicle undergo propagation -----
                rad.zmer += 1
                updateRadicle(par,rad)
                tnextP[j] += propTime
            end
            curL = multiLbounceStepUpdate(par,rad,tempτ)
            minL[i,j] = min( minL[i,j] ,curL)
            # ------  ----- ----- -----  -----
        end
    end
    return(P=pos,V=vel,Zmer=zmers,Time=Time,minL=minL)
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
