using ProgressBars
using Plots

D_mon = 1e-9; # 这个是什么？df 吗？
parDia = 90;  # Particle Diameter

maxTimeFac = (1^2)/2;   # 这又是啥？Maximum Displacement in 1-D Factor
t = 0;

dtCase1 = maxTimeFac/(D_mon*(1e14));   # Maximum jump time for the monomer
dt = min(0.01,dtCase1);  # Time taken in each step
distL = (2*D_mon*(1e14)*dt)^0.5; # 这个是什么？

numSimuPoint = Int(ceil(1/dt));
storeLoc = zeros(numSimuPoint+1,5);

# Starting Location
theta = 0.5 * 2 * pi;  # Same Starting point in each simulation
phi = acos((2*0.5)-1);
storeLoc[1,1] = (parDia/2)* cos(theta) * sin(phi);  # X-Coordinate
storeLoc[1,2] = (parDia/2)* sin(theta) * sin(phi);  # Y-Coordinate
storeLoc[1,3] = (parDia/2)* cos(phi);
storeLoc[1,4] = (parDia/2);

for ii in ProgressBar(1:numSimuPoint)
    flagX = 0;

    while (flagX == 0)
        newLoc = storeLoc[ii,1:3] .+  (distL .* randn(3));
        newR = (sum(newLoc.^2))^0.5;
        if (newR < (parDia/2))
            flagX = 1;
            storeLoc[ii+1,1:3] = newLoc;
            storeLoc[ii+1,4] = newR;
        end
    end
    
    t = t + dt;
    storeLoc[ii+1,5] = t;
end

XXXX = zeros(numSimuPoint+1,1);
for ii = 1:3
    XXXX = XXXX + ((storeLoc[:,ii] .- storeLoc[1,ii]).^2);
end
XXXX = XXXX.^0.5;

plot(storeLoc[:,5],XXXX)
