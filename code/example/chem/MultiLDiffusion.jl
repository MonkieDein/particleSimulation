include("../../chem.jl")

# Diffusion coefficient parameters
Wp = 1-Wm                               # Wp : Weightage of polymer
T = 80;                                 # T : Reaction Temparature (°C)
Tg_sys = 25;                            # Tg : Glass Transition Temperature (ignore plasticization w/ Monomer)
D = 10^logD(Wp ,T - Tg_sys,unit="nm")   # D : diffusion constant / coefficient
# zmer length and number of propagation 
lZmerl = 4;                             # Zmer length (Different for each monomer): 
lPropl = 14;  


