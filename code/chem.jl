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
    L :: DataFrame # Layers rstart and Diffusion coefficient
    Wp :: Float64 # Weightage of polymer
    function mLparticle(x,y,z,R,L,Wp)
        new(Round(x,y,z,R),L,Wp)
    end
    function mLparticle(P,R,L,Wp)
        new(coord(P),R,L,Wp)
    end
end

mutable struct Radicle
    obj :: Round # radicle itself as a Round object
    D :: Float64 # diffusion coefficient
    zmer :: Int # zipper mer length
    τ :: Float64 # Global Time step
    σx :: Float64 # Mean square displacement of each axis
    l :: Int # layer index
    function Radicle(x,y,z,D;r=0,vx=0,vy=0,vz=0,σx = 1,l = 1,zmer = 1,τ = 0.01)
        new(Round([x,y,z],r,V=[vx,vy,vz]),D,zmer,τ,σx,l)
    end
    function Radicle(P,D;r=0,V=[0.0,0,0],σx = 1,l = 1,zmer = 1,τ = 0.01)
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
    D = 10^logD(particle.Wp ,particle.L[radicle.l],unit="nm")
    radicle.D = D/(radicle.zmer^(0.5+1.75*particle.Wp))
    newσx = sqrt( 2 * radicle.D / radicle.τ)
    radicle.v = coord( vec(radicle.v) .* (newσx / radicle.σx) ) # scale v multiplying (new_σx/old_σx) 
    radicle.σx = newσx
end

function multiLbounceStepUpdate(particle::mLparticle,radicle::Radicle,Δt)  
    tempΔt = Δt
    while tempΔt > 0
        r = L2Distance(vec(particle.obj.p),vec(radicle.obj.p))
        I = [intersections(particle.obj, R, radicle.obj, tempΔt) for R in particle.L.R]
        if minimum(I) < 1       # Transition over Layers
            (tmin,l) = findmin(I)
            if tmin == 0        # Update parameter before transition
                d = normalize( vec(radicle.obj.p) .- vec(particle.obj.p) )
                v = normalize( vec(radicle.obj.v) ) 
                radicle.l = l - (dot(d,v) > 0) # If dot product > 0 means going toward the center
                updateRadicle(particle,radicle) # Update radicle new layer information
                tmin = I[partialsortperm(I, 2)] # find the second smallest time
            end
            updateMotion(radicle.obj,tmin)
            tempΔt -= tmin
        else                    # Collide on particle wall
            colission = collisionNreflection(particle.obj,radicle.obj,tempΔt)
            updateMotion(radicle.obj,colission.time * tempΔt)
            if (colission.time < 1)
                radicle.obj.v = colission.reflection
            end
            tempΔt -= colission.time * tempΔt
        end
    end
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



