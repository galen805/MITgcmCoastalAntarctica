% This is a matlab script that generates the input data for curtain sim

% created by gwilcox 6/21/2025

close all

%% initialize grid and time discretization
northen_end_lat = -76; % starting at -76 lat

% CENTERED Dimensions of grid, high resolution region
% uniform resolution for this model
nx=288; % total x points
nx1=nx; % central study region
ny=240;
nz=32;

hres = 2; % set resolution in km... Arthun 2013 says < 2km is required
dlat = hres/111;
dlon = hres/(111*cosd(northen_end_lat));

delta_z = 25;

acc = 'real*8';
ieee = 'b';

% Grid file
dxfile = dlon*ones(1,nx1);
dyfile = dlat*ones(1,ny);
dzfile = delta_z*ones(1,nz);

% lengthen zonal direction
% increase=0.05;
% nx_west = (nx-nx1)/2;
% nx_east = nx_west;
% 
% east_add(1)=dlon*(1+increase);
% for i=2:nx_east
%     east_add(i) = east_add(i-1)*(1+increase);
% end
% west_add = flip(east_add);
% dxfile = [west_add dxfile east_add];

% print out the bin files
fid=fopen('dxfile.bin','w','b'); fwrite(fid,dxfile,acc);fclose(fid);
fid=fopen('dyfile.bin','w','b'); fwrite(fid,dyfile,acc);fclose(fid);
fid=fopen('dzfile.bin','w','b'); fwrite(fid,dzfile,acc);fclose(fid);


%% longtitude and latitude vectors at face g and center c

latg = northen_end_lat+cumsum(dyfile);
latc = latg+dyfile/2;

long = cumsum(dxfile);
lonc = long+dxfile/2;

zgp1 = [0,cumsum(dzfile)];
zg = zgp1(1:end-1);
zc = .5*(zgp1(1:end-1)+zgp1(2:end));


%% Bathymetry parameters

% construct the center high resolution part first
% Point out the land part, bathy=0
% ice shelf part is about 50 points in x, 100 points in y

% H- Nominal depth of model (meters)
H = -500;
bathy_combined = ones(nx,ny)*H;
iceshelf_width_cells=60; % even number of grid cells
west_land_index=nx/2-iceshelf_width_cells/2;
east_land_index = nx/2+iceshelf_width_cells/2;

west_iceshelf_index = west_land_index+1;
east_iceshelf_index = east_land_index-1;

% bathy_high_reso(1,:) = 0; % Wall on the west side
bathy_combined(:,1) = 0; % Wall on the south side

iceshelf_north_end_index = 60; % important, indicate land end
%
bathy_combined(1:(west_land_index), 1:iceshelf_north_end_index) = 0;
bathy_combined((east_land_index):nx, 1:iceshelf_north_end_index) = 0;
%

%% Round up the land corners

% Point land northwest
land_radius = 15;
i_center = west_land_index-land_radius;
j_center = iceshelf_north_end_index-land_radius;
for i = i_center:(west_land_index)
    for j = j_center:iceshelf_north_end_index
        if ((i-i_center)^2+(j-j_center)^2) > land_radius^2
            bathy_combined(i,j) = H;
        end
    end
end

% Point land northeast
i_center = east_land_index+land_radius;
j_center = iceshelf_north_end_index-land_radius;
for i = (east_land_index):i_center
    for j = j_center:iceshelf_north_end_index
        if ((i-i_center)^2+(j-j_center)^2) > land_radius^2
            bathy_combined(i,j) = H;
        end
    end
end


% Point iceshelf southwest
i_center = west_iceshelf_index+land_radius;
j_center = 2+land_radius;
for i = west_iceshelf_index:i_center
    for j = 2:j_center
        if ((i-i_center)^2+(j-j_center)^2) > land_radius^2
            bathy_combined(i,j) = 0;
        end
    end
end


% Point iceshelf southeast
i_center = east_iceshelf_index-land_radius;
j_center = 2+land_radius;
for i = i_center:east_iceshelf_index
    for j = 2:j_center
        if ((i-i_center)^2+(j-j_center)^2) > land_radius^2
            bathy_combined(i,j) = 0;
        end
    end
end



%% continental shelf shape gwilcox 2/28/24

% main geometry points
% flip sign of max term to control prograde vs retrograde
% for now, prograde (-max)
Hc=500; % depth at back of cavity, m
hf=500; % depth at N edge of domain, m

y_bathy_ambient=latc-latc(1);

% tanh shape shelfbreak
shelfbreak_ind=ny;
L_shelfbreak=0.2;
ha=(Hc+hf)/2+(hf-Hc)/2*tanh((latc-latc(shelfbreak_ind))/L_shelfbreak);
ha=-ha;

% modify bathy_combined
bathy_sloped=zeros(size(bathy_combined));
for i=1:ny
    bathy_sloped(:,i)=ha(i);
end

% add to bathy_combined
bathy_sloped(bathy_combined==0)=0;
bathy_combined=bathy_sloped;


%% add trough to ambient bathymetry
% retrograde slope

% trough shape parameters
Hc_t=800; % depth at south end of domain
hf_t=650; % trough depth at north end of domain

% calculate center depth of trough

% retrograde
Ht=Hc_t+(hf_t-Hc_t)/y_bathy_ambient(shelfbreak_ind)*y_bathy_ambient;
Ht=-Ht;

% modify to match shelfbreak
Ht(Ht>ha)=ha(Ht>ha);

% smooth a few times like icetopo
% need to smooth both so they match at the offshore end!
smoothen_times = 300;
for i = 1:smoothen_times
    Ht = shapiro1(Ht, 2, 1);
    ha = shapiro1(ha, 2, 1);
end

% trough width in cells
width_cavity=iceshelf_width_cells; % in cavity, match to ice shelf width
%width_coast=30; % beyond cavity, turns into a constant width trough
width_coast=width_cavity;
width_slope=width_coast;

% convert to degrees
width_cavity_deg=dlon*width_cavity;
width_coast_deg=dlon*width_coast;
width_slope_deg=dlon*width_slope;

% set up tanh trough width function
widthchange=1.125;
widthchange_scale=.125;
y_bathy_ambient=latc-latc(1);
trough_width_deg=width_coast_deg-(width_cavity_deg-width_coast_deg)/2*(tanh((y_bathy_ambient-widthchange)/widthchange_scale)-1);

%plot(trough_width_deg)

%trough_width_deg=width_coast+(width_slope-width_coast)/y_bathy_ambient(end)*y_bathy_ambient;
Wc=trough_width_deg/4;
trough_center=(lonc(nx/2+1)+lonc(nx/2))/2*ones(size(latc));
%trough_center=((lonc(end)-lonc(1))+lonc(1))/2*ones(size(latc));

ht_depth=Ht-ha;
% add trough to bathy_combined
bathy_trough=zeros(size(bathy_combined));
for i=1:ny
    % gaussian
    % bathy_trough(:,i)=ha(i)+ht_depth(i)*exp(-(lonc-trough_center(i)).^2/Wc(i)^2);

    % steep-sided
    side_scale=0.6; % deg lat
    bathy_trough(1:nx/2,i)=ha(i)+ht_depth(i)/2*(tanh((lonc(1:nx/2)-trough_center(i)+trough_width_deg(i)/2)/side_scale)+1);
    bathy_trough(nx/2+1:end,i)=ha(i)-ht_depth(i)/2*(tanh((lonc(nx/2+1:end)-trough_center(i)-trough_width_deg(i)/2)/side_scale)-1);
end
bathy_trough(bathy_combined==0)=0; % set coastal part to zero
bathy_combined=bathy_trough;

if local_plot

    % bathy lines
    figure;
    hold on
    plot(latc,Ht,'-','linewidth',2)
    plot(latc,ha,'-','linewidth',2)
    legend('trough','nominal bathy')
    xlabel('deg lat')
    ylabel('m')
    hold off
    grid on
    xlim([latc(1),latc(end)])


    % plot trough shape
    figure;
    hold on
    % at edge of ice shelf
    shape_index=iceshelf_north_end_index+1;
    lonc_centered=lonc-trough_center(shape_index);
    xcoord_km=lonc_centered*111.32*cosd(latc(shape_index));
    plot(xcoord_km,bathy_combined(:,shape_index),'-','linewidth',2)
    yline(ha(shape_index),'--')
    xlabel('km')
    ylabel('m')  
    xlim([-100,100])
    grid on
    set(gca,'fontsize',12)
    title('Trough shape')

    figure;
    hold on

    X=lonc;
    Y=latc;
    [latc_grid,lonc_grid]=meshgrid(Y,X);

    surf(latc_grid,lonc_grid,bathy_combined)
    colorbar
    view(45, 30)

    grid on
    box on
    axis tight

end





%% Ice shelf

% A more smooth ice line with Shapiro filter

H1=-800;
tp12=10; % turning point 1
H2=-500;
tp23=iceshelf_north_end_index-5; % turning point 2
H3=-200;
ice_line=zeros(1,ny);

% grounding line
ice_line(1:tp12)=H1+(H2-H1)*(latc(1:tp12)-latg(1))/(latc(tp12)-latg(1));

% main part
ice_line(tp12:tp23)=H2+(H3-H2)*(latc(tp12:tp23)-latc(tp12))/(latc(tp23)-latc(tp12));

% pointy part
ice_line(tp23:iceshelf_north_end_index)=H3+(0-H3)*(latc(tp23:iceshelf_north_end_index)-latc(tp23))/(latc(iceshelf_north_end_index)-latc(tp23));

% Smoothen iceshelf topo for 100 times
smoothen_times = 50;
for i = 1:smoothen_times
    ice_line(1:(iceshelf_north_end_index)) = shapiro1(ice_line(1:(iceshelf_north_end_index)), 2, 1);
end

ice_line(:,(iceshelf_north_end_index+1):end) = 0; % ocean part, no iceshelf

icetopo = ones(size(bathy_combined,1),1)*ice_line;

% eliminate anywhere icetopo is deeper than bathy
% 1/23/24
bathy_combined(abs(icetopo)>abs(bathy_combined))=0;
icetopo(bathy_combined==0) = 0;

% make sure icetopo and bathy agree
shouldBe1=min(bathy_combined<=icetopo,[],'all'); % should be 1 (always deeper than icetopo)
shouldBe0=max(bathy_combined==0 & icetopo~=0,[],'all'); % should be 0 (no instances where bathy is zero (land) and ice is nonzero)

% write bathy
fid=fopen('bathy.box','w','b'); fwrite(fid,bathy_combined,acc);fclose(fid);

% write icetopo
fid=fopen('icetopo.exp1','w','b'); fwrite(fid,icetopo,acc);fclose(fid);


%% wind and air external forcings
vwind_main = [0,0]; % not using this anymore

% close all
% make katabatic wind pattern
x0=nx/2;
y0=iceshelf_north_end_index+1;
r2 = 15;   % decay scale in x
r3 = 50;   % decay scale in y
r1 = iceshelf_width_cells-r2/2;   % width of tanh region in x

vwindkat = zeros(nx, ny);
for i = 1:nx
    for j = 1:ny
        % tanh shape in x centered at x0
        tanh_x = 0.5 * (1 - tanh((abs(i - x0) - r1) / r2));
        
        % exponential decay northward from y0
        exp_y = exp(-(j - y0)/r3);
        
        vwindkat(i,j) = tanh_x * exp_y;
    end
end

% generate wind with synoptic and seasonal variability for one year
% forcing updated every 6 hours
numEntries=4*360;
vwind=zeros(nx,ny,numEntries+1);
uwind=zeros(nx,ny,numEntries+1);

% synoptic oscillation period 2 days
T=2*4;
time_6hr=0:numEntries;
frac_var=0.5; % variability as fraction of mean
synoptic_oscillation=(1-frac_var*cos(2*pi*time_6hr/T));

for t=1:numEntries+1
    % add annual oscillation between winter and summer mean values period 365 days
    vwind_kat_max(t) = mean(vwind_kat) + 0.5*(vwind_kat(1)-vwind_kat(2))*sin(2*pi*(t-1)/(360*4));
    vwind_rest(t) = mean(vwind_main) + 0.5*(vwind_main(1)-vwind_main(2))*sin(2*pi*(t-1)/(360*4));
    uwind_rest(t) = mean(uwind_main) + 0.5*(uwind_main(1)-uwind_main(2))*sin(2*pi*(t-1)/(360*4));

    % add synoptic var
    vwind_kat_max(t) = vwind_kat_max(t)*synoptic_oscillation(t);
    vwind_rest(t) = vwind_rest(t)*synoptic_oscillation(t);
    uwind_rest(t) = uwind_rest(t)*synoptic_oscillation(t);

    % combine
    uwind(:,:,t)=uwind_rest(t);
    vwind(:,:,t) = vwindkat.*vwind_kat_max(t) + vwind_rest(t);
end

% take off land part
uwind(:,1:iceshelf_north_end_index,:)=0;
vwind(:,1:iceshelf_north_end_index,:)=0;

% write bin files
fid=fopen('vwind.bin','w','b'); fwrite(fid,vwind,acc);fclose(fid);
fid=fopen('uwind.bin','w','b'); fwrite(fid,uwind,acc);fclose(fid);

%% Airtemp, LW radiation, SW radiation, humidity updated every 15 days
numATEntries=360/15;
lwdownfield_main = [200,250]; % winter, summer
swdownfield_main = [0,250]; % winter, summer
humidfield_value = 0.0005; % constant

airtemp=ones(nx,ny,numATEntries+1);
lwdownfield=ones(nx,ny,numATEntries+1);
swdownfield=ones(nx,ny,numATEntries+1);
humidfield=ones(nx,ny,numATEntries+1);
for t=1:numATEntries+1
    % annual trends
    airtemp_val(t)=mean(airtemp_main) + 0.5*(airtemp_main(1)-airtemp_main(2))*sin(2*pi*(t-1)/(360/15)); 
    airtemp(:,:,t)=airtemp_val(t);

    lwdownfield_val(t)=mean(lwdownfield_main) + 0.5*(lwdownfield_main(1)-lwdownfield_main(2))*sin(2*pi*(t-1)/(360/15));
    lwdownfield(:,:,t)=lwdownfield_val(t);

    swdownfield_val(t)=mean(swdownfield_main) + 0.5*(swdownfield_main(1)-swdownfield_main(2))*sin(2*pi*(t-1)/(360/15));
    swdownfield(:,:,t)=swdownfield_val(t);

    humidfield(:,:,t)=humidfield_value;

end

fid=fopen('airtempfield.bin','w','b'); fwrite(fid,airtemp,acc);fclose(fid);
fid=fopen('lwdownfield.bin','w','b'); fwrite(fid,lwdownfield,acc);fclose(fid);
fid=fopen('swdownfield.bin','w','b'); fwrite(fid,swdownfield,acc);fclose(fid);
fid=fopen('humidfield.bin','w','b'); fwrite(fid,humidfield,acc);fclose(fid);

%% Sea ice initial files
SIarea_init_value = 1;
SIarea_init = SIarea_init_value*ones(size(bathy_combined));
SIarea_init(bathy_combined == 0) = 0; % get rid of the land
SIarea_init(icetopo ~= 0) = 0; % get rid of ice shelf, ice tongue and pack ice
fid=fopen('SIarea_init.bin','w','b'); fwrite(fid,SIarea_init,acc);fclose(fid);

SIheff_init_value = 0.5; % meter
SIheff_init = SIheff_init_value*ones(size(bathy_combined));
SIheff_init(bathy_combined == 0) = 0; % get rid of the land
SIheff_init(icetopo ~= 0) = 0; % get rid of ice shelf, ice tongue and pack ice
fid=fopen('SIheff_init.bin','w','b'); fwrite(fid,SIheff_init,acc);fclose(fid);

%% Empty external forcing- use 0 eyerywhere
empty2dfield1 = 0.0*ones(size(bathy_combined));
empty2dfield2 = empty2dfield1;

empty2dfield(:,:,1) = empty2dfield1;
empty2dfield(:,:,2) = empty2dfield2;

fid=fopen('empty2dfield.bin','w','b'); fwrite(fid,empty2dfield,acc);fclose(fid);


%% sea ice zero-thickness OBCS at N boundary only

iceOBCSanorth1 = 0.0*ones(size(bathy_combined,1),1);
iceOBCSanorth2 = iceOBCSanorth1;
iceOBCSanorth(:,1) = iceOBCSanorth1;
iceOBCSanorth(:,2) = iceOBCSanorth2;
fid=fopen('iceOBCSanorth.bin','w','b'); fwrite(fid,iceOBCSanorth,acc);fclose(fid);

iceOBCShnorth = 0.0*iceOBCSanorth;
fid=fopen('iceOBCShnorth.bin','w','b'); fwrite(fid,iceOBCShnorth,acc);fclose(fid);

% east and west obcs for ice
% iceOBCSawest1 = 0.0*ones(size(bathy_combined,2),1);
% iceOBCSawest2 = iceOBCSawest1;
% iceOBCSawest(:,1) = iceOBCSawest1;
% iceOBCSawest(:,2) = iceOBCSawest2;
% fid=fopen('iceOBCSawest.bin','w','b'); fwrite(fid,iceOBCSawest,acc);fclose(fid);
% 
% iceOBCShwest = 0.0*iceOBCSawest;
% fid=fopen('iceOBCShwest.bin','w','b'); fwrite(fid,iceOBCShwest,acc);fclose(fid);

% velocity?
% uiceOBCSwest = 0.0*iceOBCSawest;
% fid=fopen('uiceOBCSwest.bin','w','b'); fwrite(fid,uiceOBCSwest,acc);fclose(fid);


%% Ptracers for meltwater and salt (tied to shelfice and seaice packages)
% Add the parts for ptracers surface intial file

ptracers2d=ones(size(bathy_combined));
ptracers2d(:,1:iceshelf_north_end_index) = 0;
ptracers3d = zeros(size(bathy_combined,1),size(bathy_combined,2),nz);

fid=fopen('ptracers_file.bin','w',ieee); fwrite(fid,ptracers3d,acc); fclose(fid);

% Add zero mask file for Ptracers, use the below for RBCS
ptracers_mask = zeros(size(ptracers3d));

fid=fopen('rbcs_mask_ptracers.bin','w',ieee); fwrite(fid,ptracers_mask,acc); fclose(fid);

%% generate 1D background and forcing ocean profiles
% CDW blend function
CDWblend = 0.5 * (1 + tanh((zc - CDW0_TD)/CDW0_dz_transition));  % 0 above, 1 below

% basic linear T, rho fields
dens_range = dens_top + (dens_bottom - dens_top)/zc(end)*zc;
T_range = T0*ones(1,nz);

forcingdens_range=dens_range;
forcingT_range=T_range;

% add CDW to initial condition?
if warm_ocean
    % add thermocline/CDW layer
    T_range=(1 - CDWblend).*T_range + CDWblend.*CDW0_temp;
end

% add CDW to forcing?
if warm_forcing
    % add thermocline/CDW layer
    forcingT_range=(1 - CDWblend).*forcingT_range + CDWblend.*CDW0_temp;
end

% calc salinity and strat
S_range=calc_salinity_native(dens_range,T_range);
N2_range = calc_N2(dens_range, zc);

forcingS_range=calc_salinity_native(forcingdens_range,forcingT_range);
forcingN2_range = calc_N2(forcingdens_range, zc);


% plot initial and forcing profiles
if local_plot
    figure;
    set(gcf,'pos',[300,300,600,400])
    tiledlayout(1,2);
    nexttile
    hold on
    plot(T_range,-zc,'-k','linewidth',2)
    plot(forcingT_range,-zc,'-r','linewidth',2)
    legend('Initial','OBCS')
    title('\Theta')
    ylim([-1*zc(end),-1*zc(1)])
    set(gca,'fontsize',18)
    grid on
    box on
    ylabel ('Depth (m)')
    xlabel('Deg C')

    nexttile
    hold on
    plot(S_range,-zc,'-k','linewidth',2)
    plot(forcingS_range,-zc,'-r','linewidth',2)
    legend('Initial','OBCS')
    title('SP')
    ylim([-1*zc(end),-1*zc(1)])
    set(gca,'fontsize',18)
    grid on
    box on
    xlabel('PSU')

    % nexttile
    % hold on
    % plot(dens_range,-zc,'-k','linewidth',2)
    % plot(forcingdens_range,-zc,'-r','linewidth',2)
    % legend('Initial','Forcing')
    % title('\rho_\theta')
    % ylim([-1*zc(end),-1*zc(1)])
    % set(gca,'fontsize',18)
    % grid on
    % box on
    % xlabel('kg/m^3')
    % 
    % nexttile
    % hold on
    % plot(N2_range,-zc,'-k','linewidth',2)
    % plot(forcingN2_range,-zc,'-r','linewidth',2)
    % legend('Initial','Forcing')
    % title('N^2')
    % ylim([-1*zc(end),-1*zc(1)])
    % xlim([0,3e-5])
    % set(gca,'fontsize',18)
    % grid on
    % box on
    % xlabel('1/s^2')

    exportgraphics(gcf,'TS_prof.png','resolution',300)
end

%% write initial conditions
% set up and write reference files
S_ref=S_range;
T_ref=T_range;
fid=fopen('salt_ref_file.bin','w','b'); fwrite(fid,S_ref,acc);fclose(fid);
fid=fopen('T_ref_file.bin','w','b'); fwrite(fid,T_ref,acc);fclose(fid);

% set up and write initial conditions
init_S_file = repmat(reshape(S_range, [1 1 nz]), nx, ny, 1);
init_T_file = repmat(reshape(T_range, [1 1 nz]), nx, ny, 1);
fid=fopen('init_S_file.bin','w','b'); fwrite(fid,init_S_file,acc);fclose(fid);
fid=fopen('init_T_file.bin','w','b'); fwrite(fid,init_T_file,acc);fclose(fid);


%% write OBCS
OBNtFile = repmat(reshape(forcingT_range, [1 nz]), nx, 1);
OBNsFile = repmat(reshape(forcingS_range, [1 nz]), nx, 1);

fid=fopen('OBCSnorth_S.bin','w','b'); fwrite(fid,OBNsFile,acc);fclose(fid);
fid=fopen('OBCSnorth_T.bin','w','b'); fwrite(fid,OBNtFile,acc);fclose(fid);



if local_plot

    % temp xz
    figure('pos',[200,200,600,400]);
    width=7; % relative to colorbar
    t1=tiledlayout(1,width);
    nexttile(t1,1,[1 width-1]);

    levels=[-2:0.25:1.5];
    numLevels=length(levels);
    cmap=cmocean('thermal',numLevels);

    hold on
    [lonc_grid_xz,Zcgrid_xz]=meshgrid(lonc,-zc);
    contourf(lonc_grid_xz, Zcgrid_xz, squeeze(OBNtFile(:,:))', levels,'LineColor','none')
   
    bathy_line = bathy_combined(bathy_combined(:,ny)~=0, ny)';
    lonc_bathy = lonc(bathy_combined(:,ny)~=0)' - 0.04;
    plot(lonc_bathy, bathy_line,'-k','linewidth',2); hold on;

    title('Temp')
    set(gca,'fontsize',20)

    add_discrete_colorbar(t1,width,levels,cmap,'Deg C',2,2)
    saveas(gcf,'obcs_T_xz.png')

    % temp xz
    figure('pos',[200,200,600,400]);
    width=7; % relative to colorbar
    t1=tiledlayout(1,width);
    nexttile(t1,1,[1 width-1]);

    levels=[33.5:0.1:35.1];
    numLevels=length(levels);
    cmap=cmocean('haline',numLevels);

    hold on
    [lonc_grid_xz,Zcgrid_xz]=meshgrid(lonc,-zc);
    contourf(lonc_grid_xz, Zcgrid_xz, squeeze(OBNsFile (:,:))', levels,'LineColor','none')
   
    bathy_line = bathy_combined(bathy_combined(:,ny)~=0, ny)';
    lonc_bathy = lonc(bathy_combined(:,ny)~=0)' - 0.04;
    plot(lonc_bathy, bathy_line,'-k','linewidth',2); hold on;

    title('Salinity')
    set(gca,'fontsize',20)

    add_discrete_colorbar(t1,width,levels,cmap,'PSU',2,2)
    saveas(gcf,'obcs_S_xz.png')
end

%% write to plotting_data.mat
if local_plot
save('plotting_data.mat','bathy_combined','icetopo','lonc','latc','vwindkat','zc','trough_width_deg','dxfile','dyfile') 
copyfile('plotting_data.mat', '/Users/galenwilcox/Desktop/walker_zhang_research/transition/MITgcm/transition_processing/plotting_data.mat')
end

%% clean bin files
if local_plot
    currentDir = pwd;
    files = dir(fullfile(currentDir, '*.bin'));
    for i = 1:length(files)
        filename = fullfile(currentDir, files(i).name);
        delete(filename);
    end
end

