
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

parDia = 100;   # Particle Diameter

phiM = 0.05;
monConc = (phiM/100)*1000;  # (phiM/Mw)*1000, SO CHANGE FOR DIFFERENT MONOMER
phiP = 1- phiM;

## Per Z condition monConc = 3.8 mol/L BMA -- DIRECT ENTERED BY CLIFF as phiM = 0.05

monRMS = 0.69;  

T = 80;         # Reaction Temparature (C)
Tg_sys = 25;   # Polymer Tg (without considering the plasticization with Monomer)

T = T + 273.15;
Tg_sys = Tg_sys + 273.15;

## Constants Associated with the simulation
R = 8.314;      # Gas Constnat (J/K-mol)
N_av = 6.023e23;    # Avagadro Number (1/mol)

# Propagation Rate Coefficient for Monomer
KpArA_chem = (2.673e6)*exp((-1)*(22.36e3)/(R*T));   # FOR MMA
#KpArA_chem = 1515;   # FOR nBMA (DIRECT FROM CLIFF) mol-1.s-1

propFreq = KpArA_chem*monConc;
propFreq = 1/propFreq;  # Propagation Frequency

sol = diffusion_limit(T - Tg_sys);
D_lim = sol.D_lim
D_coeff = sol.D_coeff

# Calculate Diffusion Coefficient
if (phiP > 0)
    diff_curr_l = searchsortedlast(D_lim,phiP)
else
    diff_curr_l = 1;
end
diff_coeff_curr = D_coeff[diff_curr_l,:];

D_mon = 10^(diff_coeff_curr[1] + diff_coeff_curr[2]*phiP);  # Monomer Diff

zmerLength = 4; ## CHANGE IF USING DIFFERENT MONOMER
#storeP = zeros(1000,5);

##

maxTimeFac = (1^2)/2;   # Maximum Displacement in 1-D Factor
t = 0;
tnextP = propFreq;

dtCase1 = maxTimeFac/(D_mon*(1e14));   # Maximum jump time for the monomer
dt = min(0.01,dtCase1);  # Time taken in each step
distL = (2*D_mon*(1e14)*dt)^0.5; 

numSimuPoint = 14; # number of propagation steps (after z-mer length)
storeLoc = zeros(numSimuPoint,8);
# TMP - AMIT
storeLoc_all = zeros(1000000,4);
tmp_ind = 1

# Starting Location
theta = 0.5 * 2 *pi();  # Same Starting point in each simulation (i.e. the outer surface of the particle)
phi = acos((2*0.5)-1);
storeLoc[1,1] = (parDia/2)* cos(theta) * sin(phi);  # X-Coordinate
storeLoc[1,2] = (parDia/2)* sin(theta) * sin(phi);  # Y-Coordinate
storeLoc[1,3] = (parDia/2)* cos(phi);
storeLoc[1,4] = (parDia/2);
storeLoc[1,6] = zmerLength;
newLocCurr = storeLoc[1,1:3];
newRCurr = parDia/2;

DmonCurr = D_mon/(zmerLength^(0.5+1.75*phiP));

storeLoc[1,7] = DmonCurr;
storeLoc[1,8] = 0; # This is just to set the starting distance at the particle periphery

# TMP - AMIT
storeLoc_all[tmp_ind,:] = storeLoc[1,1:4];
tmp_ind = tmp_ind + 1;

#currCount = 1;
ii = 1;

while(ii < numSimuPoint)
    
    dtCase1 = maxTimeFac/(DmonCurr*(1e14)); # Max jump time ########WHY IS THIS ALWAYS 0.5??
    dt = min(propFreq,dtCase1); 
    distL = (2*DmonCurr*(1e14)*dt)^0.5; ### IS THIS THE DISTANCE TRAVELED IN EACH PROPAGATION STEP?? WHY 2 AND NOT 6?? ###
    ######## WHY IS distL ALWAYS 1.0 WHEN THE EQUATION ABOVE DOES NOT CALCULATE TO 1.0
    ######## SHOULD distL be defined as #1.4f instead of #d or something???
    
    flagX = 0;
    while (flagX == 0)
        newLoc = newLocCurr +  (distL .* randn(3));
        newR = (sum(newLoc.^2))^0.5;
        
        if (newR < (parDia/2))
            flagX = 1;
        end
    end
    
    newLocCurr = newLoc;
    newRCurr = newR;
	
	# TMP - AMIT
	storeLoc_all[tmp_ind, 1:3] = newLoc;
	storeLoc_all[tmp_ind, 4] = newR;
	tmp_ind = tmp_ind + 1;
    
    t = t + dt;
    if (t >= tnextP)
        tnextP = tnextP + propFreq;
        ii = ii + 1;
        
        storeLoc[ii,1:3] = newLoc;
        storeLoc[ii,4] = newR;
        storeLoc[ii,5] = t;
        storeLoc[ii,6] = storeLoc[ii-1,6] + 1;
        DmonCurr = D_mon/(storeLoc[ii,6]^(0.5+1.75*phiP));
        storeLoc[ii,7] = DmonCurr;
        storeLoc[ii,8] = storeLoc[ii-1,8] + (abs(storeLoc[ii,4] - storeLoc[ii-1,4])); ## TRYING TO WRITE THE DISTANCE TRAVELED IN EACH STEP
    end
end

XXXX = zeros(numSimuPoint,1);
for ii = 1:3
    XXXX = XXXX + ((storeLoc[:,ii]-storeLoc[1,ii]).^2);
end
XXXX = XXXX.^0.5;

storeLoc_all = storeLoc_all[1:tmp_ind-1,:];

#hold on;
#plot(storeLoc(:,5),XXXX,'g:');
#save(['DiffuseProp_MMA_Tg_' num2str(Tg_sys-273.15) '_' datestr(now,'mm-dd-YY_HH-MM') '.mat'])
