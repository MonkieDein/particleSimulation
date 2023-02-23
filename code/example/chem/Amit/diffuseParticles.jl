using Random
using Plots
using Statistics

function diffusion_limit(delT)
    A1 = -4.428 + delT*0 # C1+delT*C3
    B1 = 1.842 - delT*(8.12e-3) # C2-delT*C4
    A2 = 26 + delT*0.0797
    B2 = 37
    A3 = 159 + delT*0.3664
    B3 = 170
    A4 = (-13.7) + delT*0
    B4 = 0.500
    D_lim = [0; (A2-A1)/(B2-B1);(A3-A2)/(B3-B2);(A4-A3)/(B4-B3)]
    D_coeff = [A1 (-B1);A2 (-B2);A3 (-B3);A4 (-B4)]
    return(D_lim=D_lim,D_coeff=D_coeff)
end

T = 70; # C
Tg_sys = 90; # C
phiP = 0.95;

sol = diffusion_limit(T - Tg_sys);
D_lim = sol.D_lim
D_coeff = sol.D_coeff

# Calculate Diffusion Coefficient
if (phiP > 0)
    diff_curr_l = searchsortedlast(D_lim,phiP);
else
    diff_curr_l = 1;
end
diff_coeff_curr = D_coeff[diff_curr_l,:];

D = 10^(diff_coeff_curr[1] + diff_coeff_curr[2]*phiP);  # Monomer Diff (cm^2/s)
D = D * (1e14);     # nm^2/s

# Simulation Parameter
dt = 0.0000001;
NParticles = 10000;
t = 0;
NSteps = 1000;
time = zeros(NSteps,1);
particle = zeros(NParticles,3);
Xsqr = zeros(NSteps,1);

for i in 1:NSteps
    particle = particle + ((2*D*dt)^0.5) .* randn(NParticles,3);
    XX = zeros(NParticles,1);
    for k = 1:3
        XX = XX + (particle[:,k].^2);
    end
    Xsqr[i] = sum(XX) / NParticles;
    t = t + dt;
    time[i] = t;
end

plot(time, Xsqr)
