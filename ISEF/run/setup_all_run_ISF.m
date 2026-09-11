clear; clc; close all; 

% flag to plot and check fields on local computer. Turn off for poseidon
local_plot=1;

% set gendata parameters
wind_main = 20;
uwind_value=[-5,-5];
linstrat=1; % controls salinity and temp stratification structure
salt_top=33.25;
%salt_bottom=34.5; % kept constant
% T_top, T_bottom set in gendata


% save wind file for gendata
% wind_scale_6hr_200days_sporadic; 
wind_scale_6hr_400days_sporadic;


% generate bin files
gendata_ISF;
shelfice_binfile;

% plot relevant model setup data
if local_plot
    plot_wind_fluc_time; 
    plot_model_setup_addUwind;
    
    paper_setup_plots6
end

