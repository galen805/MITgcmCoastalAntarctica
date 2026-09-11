clear; clc; close all; 

% flag to plot and check fields on local computer. Turn off for poseidon
local_plot=1;

% atmospheric parameters 
% entered as [winter, summer]
vwind_kat = [30,10]
uwind_main = [-10,-5]
airtemp_main = [-25,-10]

% generate bin files
gendata;
shelfice_binfile;

% plot relevant model setup data
if local_plot  
    model_setup_v2
end
