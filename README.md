# PLOSBIO_Granato_E_EtAl

All modeling data and code. Code plot outputs additionally formatted in Adobe Illustrator for publication. 
Download matlab script files into a single folder for ideal run conditions. 

All scripts require the function script HGT_func_3.m (ODE's)
They also require HGT_ss_3.m (Steady state determination)

Figures 1, 3, and S5 all use the same core scripts. There are 3 of them. 
They are: 
Fig_abc.m
Fig_def.m
Fig_ghi.m

There are initial boolean logic values to set to 1 or 0 depending which plot you wish to generate. It is all clearly annotated in the code. For example, in the Fig_abc.m script, it reads:

%To generate Figure 1: Private_Attacker = 0 Private_Target = 0;
%To generate Figure 3: Private_Attacker = 1 Private_Target = 1;
%To generate Figure S5: Private_Attacker = 0 Private_Target = 1;    

% Set 1 = True. 0 = False. If both = 1 Then Fig_c
conjugation =   1; %if 1 and toxins = 0. Then Fig_b
toxins =        0; %if 1 and conjugation = 0. Then Fig_a

% Set 1 = True. 0 = False
%Private Nutrients. Both 0 = Fig 1.
Private_Attacker = 1;
Private_Target = 1;

Figure S8 = Parameter_sweeps_density_x_frequency.m
Figure S9 = All other Parameter_sweeps_[].m files

Modeling_RawData contains .csv of all simulation data. 
