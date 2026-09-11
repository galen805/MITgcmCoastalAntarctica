clear; clc; %close all; 

% flag to plot and check fields on local computer. Turn off for poseidon
local_plot=1;

% experiment flags
warm_ocean = 0 % CDW present in initial condition
warm_forcing = 1 % CDW present in N boundary OBCS

warm_wind = 0 % warm atmospheric forcing

%% initial/boundary conditions
% basic properties
T0 = -1.8;
% dens_top = 27.5;
dens_top = 27.2;
% dens_top = 26.9;
dens_bottom = 27.8;

% CDW properties (set OFF the shelf, shouldnt vary)
CDW0_temp = 1.5;
CDW0_TD = 400;
CDW0_dz_transition = 50;

%% atmospheric forcing
% entered as [winter, summer]
uwind_main = [-5,0]; % controls heat flux, held constant

if warm_wind % warm shelf atmo
    vwind_kat = [10,0] % added directly on top of main vwind
    airtemp_main = [-10,0]
else % dense shelf atmo
    vwind_kat = [25,5] % added directly on top of main vwind
    airtemp_main = [-35,-5]
end

%% generate bin files
gendata;
shelfice_binfile;

% plot relevant model setup data
if local_plot  
    model_setup_v3
end
