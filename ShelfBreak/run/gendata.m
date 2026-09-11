% This is a matlab script that generates the input data for curtain sim

% created by gwilcox 6/21/2025

close all

%% initialize grid and time discretization
northen_end_lat = -76; % starting at -76 lat

% CENTERED Dimensions of grid, high resolution region
% uniform resolution for this model
nx=396;
ny=288;
nz=80;


hres = 2; % can set resolution in km if you prefer... Arthun 2013 says < 2km is required
dlat = hres/111;
dlon = hres/(111*cosd(northen_end_lat));

% dlat = 0.01; % manual approach
% dlon = 0.04;
delta_z = 25;

acc = 'real*8';
ieee = 'b';

% time step
dt=60;

% Grid file
dxfile = dlon*ones(1,nx);
dyfile = dlat*ones(1,ny);
dzfile = delta_z*ones(1,nz);

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
hf=2000; % depth at N edge of domain, m

y_bathy_ambient=latc-latc(1);

% tanh shape shelfbreak
shelfbreak_ind=200;
L_shelfbreak=0.15;
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
hf_t=600; % trough depth at north end of domain

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
    saveas(gcf,'trough_section.png')

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

% Smoothen iceshelf topo for 50 times
smoothen_times = 50;
for i = 1:smoothen_times
    ice_line(1:(iceshelf_north_end_index-3)) = shapiro1(ice_line(1:(iceshelf_north_end_index-3)), 2, 1);
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
% northwall_bathy = bathy_combined;
% northwall_bathy(:,end) = 0;

fid=fopen('bathy.box','w','b'); fwrite(fid,bathy_combined,acc);fclose(fid);

% write icetopo
fid=fopen('icetopo.exp1','w','b'); fwrite(fid,icetopo,acc);fclose(fid);


%% wind and air external forcings
vwind_main = [0,0]; % not using this anymore

% close all
% make katabatic wind pattern
x0=nx/2;
y0=iceshelf_north_end_index+1;
r1 = iceshelf_width_cells;   % width of tanh region in x
r2 = 20;   % decay scale in x
r3 = 80;   % decay scale in y

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

% make uwind WSC pattern
WSC_width = 50; % cells
WSC_ind = ny-WSC_width;
decayto_percent = 0; % eg, -1 is decay from uwindmain to -uwindmain
tanh_wsc = 0.5 * (1 - tanh(([1:ny]-WSC_ind)/WSC_width));  % 1 at WSC_ind, 0 far north
tanh_wsc = decayto_percent + (1 - decayto_percent) * tanh_wsc;
uwindpattern=repmat(reshape(tanh_wsc,[1,ny]),nx,1);

uwindpattern=ones(nx,ny); % or undo that if you want uniform wind



% generate wind with synoptic and seasonal variability for one year
% forcing updated every 6 hours, 1+4*360=1441
% uwind=repmat(uwind1,1,1,1441);
vwind=zeros(nx,ny,1441);
uwind=zeros(nx,ny,1441);

% synoptic oscillation period 2 days
T=2*4;
time_6hr=0:1440;
frac_var=0.5; % variability as fraction of mean
synoptic_oscillation=(1-frac_var*cos(2*pi*time_6hr/T));

for t=1:1441
    % annual oscillation between winter and summer mean values period 365 days
    vwind_kat_max(t) = mean(vwind_kat) + 0.5*(vwind_kat(1)-vwind_kat(2))*sin(2*pi*(t-1)/(360*4));
    vwind_rest(t) = mean(vwind_main) + 0.5*(vwind_main(1)-vwind_main(2))*sin(2*pi*(t-1)/(360*4));
    uwind_rest(t) = mean(uwind_main) + 0.5*(uwind_main(1)-uwind_main(2))*sin(2*pi*(t-1)/(360*4));
    
    % add synoptic var
    vwind_kat_max(t) = vwind_kat_max(t)*synoptic_oscillation(t);
    vwind_rest(t) = vwind_rest(t)*synoptic_oscillation(t);
    uwind_rest(t) = uwind_rest(t)*synoptic_oscillation(t);

    % combine
    uwind(:,:,t)=uwindpattern.*uwind_rest(t);
    vwind(:,:,t) = vwindkat.*vwind_kat_max(t) + vwind_rest(t);
end

% take off land part
uwind(:,1:iceshelf_north_end_index,:)=0;
vwind(:,1:iceshelf_north_end_index,:)=0;

% write bin files
fid=fopen('vwind.bin','w','b'); fwrite(fid,vwind,acc);fclose(fid);
fid=fopen('uwind.bin','w','b'); fwrite(fid,uwind,acc);fclose(fid);

%% Airtemp, LW radiation, SW radiation, humidity updated every 15 days
lwdownfield_main = [200,250]; % winter, summer
swdownfield_main = [0,250]; % winter, summer
humidfield_value = 0.0005; % constant

airtemp=ones(nx,ny,25);
lwdownfield=ones(nx,ny,25);
swdownfield=ones(nx,ny,25);
humidfield=ones(nx,ny,25);
for t=1:25
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

SIheff_init_value = 1.0; % meter
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


%% sea ice zero-thickness OBCS files
% N boundary
iceOBCSanorth1 = 0.0*ones(size(bathy_combined,1),1);
iceOBCSanorth2 = iceOBCSanorth1;
iceOBCSanorth(:,1) = iceOBCSanorth1;
iceOBCSanorth(:,2) = iceOBCSanorth2;
fid=fopen('iceOBCSanorth.bin','w','b'); fwrite(fid,iceOBCSanorth,acc);fclose(fid);

iceOBCShnorth = 0.0*iceOBCSanorth;
fid=fopen('iceOBCShnorth.bin','w','b'); fwrite(fid,iceOBCShnorth,acc);fclose(fid);

% % E and W boundaries
% iceOBCSawest1 = 0.0*ones(size(bathy_combined,2),1);
% iceOBCSawest2 = iceOBCSawest1;
% iceOBCSawest(:,1) = iceOBCSawest1;
% iceOBCSawest(:,2) = iceOBCSawest2;
% fid=fopen('iceOBCSawest.bin','w','b'); fwrite(fid,iceOBCSawest,acc);fclose(fid);
% 
% iceOBCShwest = 0.0*iceOBCSawest;
% fid=fopen('iceOBCShwest.bin','w','b'); fwrite(fid,iceOBCShwest,acc);fclose(fid);
% 
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

%% construct offshore temp and density profiles (ACC conditions)
close all
% generate temp profile
T_top=-1.6;
T_TD1=-1.6;
T_TD2=1.5;
T_bottom=0.5;

TD = 400;
dz_transition = 100;  % half-width of transition

% construct and blend profiles
upper_range=T_top+(T_TD1-T_top)/TD*zc;
lower_range=T_TD2+(T_bottom-T_TD2)/(zc(end)-TD)*(zc-TD); % lower linear part
blend = 0.5 * (1 + tanh((zc - TD)/dz_transition));  % 0 above, 1 below
T_range=(1 - blend).*upper_range + blend.*lower_range;
T_range_offshore=T_range;

% generate density profile
dens_top=27.5;
dens_TD1=27.6;
dens_TD2=27.8;
dens_bottom=28.0;

PC = 400;
dz_transition = 100;  % half-width of transition

% construct and blend profiles
upper_range=dens_top+(dens_TD1-dens_top)/PC*zc;
lower_range=dens_TD2+(dens_bottom-dens_TD2)/(zc(end)-PC)*(zc-PC); % lower linear part
blend = 0.5 * (1 + tanh((zc - PC)/dz_transition));  % 0 above, 1 below
dens_range=(1 - blend).*upper_range + blend.*lower_range;
dens_range_offshore=dens_range;


% calculate strat
gravity=9.81;
rhoConst = 1030;
for k=1:length(zc)-1
    dz=-(zc(k)-zc(k+1));
    drhodz=(dens_range(k)-dens_range(k+1))/dz;
    N2_offshore(k)=-gravity/rhoConst*drhodz;
end


%% construct onshore profiles (cold ocean)
% generate temp profile
T_top=-1.6;
T_bottom=-1.9;
alpha=6;
T_range = T_bottom + (T_top-T_bottom)*exp(-alpha*zc/zc(end));
T_range_onshore=T_range;

% generate density profile
dens_top=27.8;
dens_bottom=28.0;
alpha=8;
dens_range = dens_bottom + (dens_top-dens_bottom)*exp(-alpha*zc/zc(end));
dens_range_onshore=dens_range;

% calculate strat
gravity=9.81;
rhoConst = 1030;
for k=1:length(zc)-1
    dz=-(zc(k)-zc(k+1));
    drhodz=(dens_range(k)-dens_range(k+1))/dz;
    N2_onshore(k)=-gravity/rhoConst*drhodz;
end


%% blend profiles in 2D
blend_location=shelfbreak_ind;
blend_width=0.25;
blendy=1/2*(1+tanh((latc-latc(blend_location))/blend_width));

blend_yz = repmat(reshape(blendy, [ny 1]), 1, nz);

temp_off2d = repmat(reshape(T_range_offshore, [1 nz]), ny, 1);
temp_on2d = repmat(reshape(T_range_onshore, [1 nz]), ny, 1);
dens_off2d = repmat(reshape(dens_range_offshore, [1 nz]), ny, 1);
dens_on2d = repmat(reshape(dens_range_onshore, [1 nz]), ny, 1);

T_2D = temp_on2d.*(1-blend_yz) + temp_off2d.*blend_yz;


dens_2D = dens_on2d.*(1-blend_yz) + dens_off2d.*blend_yz;


%% calculate salinity in 2D
% match that densitythroughout to create salinity field

p=0.1; % surface pressure
for l=1:ny
    for k=1:nz
        n_iter=0;
        rescheck=1;

        match_dens=dens_2D(l,k) + 1000;

        % Initial guess
        S = 33.0;

        while rescheck > 1e-5
            % Compute current density and residual
            dens = densjmd95(S, T_2D(l,k), p);
            res = dens - match_dens;

            % Finite difference approximation of d(dens)/dS
            dS = 1e-4;
            dens_dS = densjmd95(S + dS, T_2D(l,k), p);
            deriv = (dens_dS - dens) / dS;

            % Prevent division by zero or non-descent
            if abs(deriv) < 1e-10
                fprintf('Small derivative; stopping iteration\n');
                break;
            end

            % Newton-Raphson step
            S_new = S - res / deriv;

            % Update for next iteration
            rescheck = abs(res);
            S = S_new;
            n_iter = n_iter + 1;

            % Optional: prevent runaway
            if n_iter > 100
                fprintf('Too many iterations at (%d,%d,%d)\n', i, j, k);
                break;
            end
        end

        ambient_S_2D(l,k) = S;

    end
end

if local_plot
    set(groot, "defaultFigureWindowStyle", "normal");

    figure;
    title('PSU')
    [lonc_grid_yz,Zc_grid_yz]=meshgrid(latc,-zc);
    contourf(lonc_grid_yz,Zc_grid_yz,ambient_S_2D',20)
    cmocean('haline')
    colorbar
    hold on
    plot(latc,Ht,'-k','linewidth',2)
    plot(latc,ha,'--k','linewidth',2)
    legend('','trough','nominal bathy','location','sw')
    set(gca,'fontsize',20)

    figure;
    title('Theta')
    [lonc_grid_yz,Zc_grid_yz]=meshgrid(latc,-zc);
    contourf(lonc_grid_yz,Zc_grid_yz,T_2D',20)
    cmocean('thermal')
    colorbar
    hold on
    plot(latc,Ht,'-k','linewidth',2)
    plot(latc,ha,'--k','linewidth',2)
    legend('','trough','nominal bathy','location','sw')
    set(gca,'fontsize',20)

    figure;
    title('rho')
    [lonc_grid_yz,Zc_grid_yz]=meshgrid(latc,-zc);
    contourf(lonc_grid_yz,Zc_grid_yz,dens_2D',20)
    cmocean('dense')
    colorbar
    hold on
    plot(latc,Ht,'-k','linewidth',2)
    plot(latc,ha,'--k','linewidth',2)
    legend('','trough','nominal bathy','location','sw')
    set(gca,'fontsize',20)
end


%% full bin files
for i=1:length(lonc)
    ambient_temp(i,:,:)=T_2D;
    ambient_S(i,:,:)=ambient_S_2D;
end
init_T_file=ambient_temp;
init_S_file=ambient_S;

% still using these for initial conditions
fid=fopen('init_S_file.bin','w','b'); fwrite(fid,init_S_file,acc);fclose(fid);
fid=fopen('init_T_file.bin','w','b'); fwrite(fid,init_T_file,acc);fclose(fid);


% write reference files
T_ref=T_range_onshore;
salt_ref=squeeze(ambient_S_2D(1,:));
fid=fopen('T_ref_file.bin','w','b'); fwrite(fid,T_ref,acc);fclose(fid);
fid=fopen('salt_ref_file.bin','w','b'); fwrite(fid,salt_ref,acc);fclose(fid);


%% create RBCS mask region for eastern side and save RBCS files
% mask specification
max_restore=1; % max restoration every N days
restore_scale=20; % half width of transition
rbcs_region_width = 30;
rbcs_y_ramp=1/2*(1+tanh(([1:ny]-(ny-rbcs_region_width))/restore_scale));
blend_rbcs_y = repmat(reshape(rbcs_y_ramp, [1 ny 1]), nx, 1, nz);


restore_top = 1;
restore_TD1 = 1;
restore_TD2 = 1; 
restore_bottom = 0;
restore_end = 1600;
upper_range=restore_top+(restore_TD1-restore_top)/TD*zc;
lower_range=restore_TD2+(restore_bottom-restore_TD2)/(restore_end-TD)*(zc-TD); % lower linear part
lower_range(lower_range<0)=0;
blend = 0.5 * (1 + tanh((zc - TD)/dz_transition));  % 0 above, 1 below
rbcs_z_ramp = (1 - blend).*upper_range + blend.*lower_range;
smoothen_times = 100;
for i = 1:smoothen_times
    rbcs_z_ramp = shapiro1(rbcs_z_ramp, 2, 1);
end

blend_rbcs_z = repmat(reshape(rbcs_z_ramp, [1 1 nz]), nx, ny, 1);

blend_rbcs = blend_rbcs_y.*blend_rbcs_z;

rbcs_mask_file = 1/24/max_restore*blend_rbcs;


if local_plot
    % plot restoration region
    figure('pos',[200,200,700,600]);
    width=7; % relative to colorbar
    t1=tiledlayout(1,width);
    nexttile(t1,1,[1 width-1]);

    levels=[0:2:20];
    numLevels=length(levels);
    cmap=flipud(cmocean('matter',numLevels));

    yzdat=squeeze(rbcs_mask_file(nx/2,:,:))';
    yzdat = 1./yzdat/24;

    [latc_grid_yz,Zc_grid_yz]=meshgrid(latc,-zc);
    contourf(latc_grid_yz,Zc_grid_yz,yzdat,levels)
    hold on
    plot(latc,Ht,'-k','linewidth',2)
    plot(latc,ha,'--k','linewidth',2)
    title('Restoration timescale')
    set(gca,'fontsize',20)
    box on

    add_discrete_colorbar(t1,width,levels,cmap,'Days',2,2)

    saveas(gcf,'restoration_region.png')
end

% write mask file
fid=fopen('rbcs_mask_file.bin','w','b'); fwrite(fid,rbcs_mask_file,acc);fclose(fid);

fid=fopen('rbcs_S_file.bin','w','b'); fwrite(fid,init_S_file,acc);fclose(fid);
fid=fopen('rbcs_T_file.bin','w','b'); fwrite(fid,init_T_file,acc);fclose(fid);


%% plot the various profiles
if local_plot
    figure;
    set(gcf,'pos',[300,300,1400,800])
    tiledlayout(2,4);
    nexttile
    hold on
    plot(T_range_offshore,-zc,'-','linewidth',2)
    title('\Theta offshore')
    set(gca,'fontsize',18)
    grid on
    box on
    ylabel ('Depth (m)')
    xlabel('Deg C')

    nexttile
    hold on
    plot(squeeze(ambient_S_2D(end,:)),-zc,'-','linewidth',2)
    title('SP offshore')
    set(gca,'fontsize',18)
    grid on
    box on
    xlabel('PSU')

    nexttile
    hold on
    plot(dens_range_offshore,-zc,'-','linewidth',2)
    title('\rho_\theta offshore')
    set(gca,'fontsize',18)
    grid on
    box on
    xlabel('kg/m^3')

    nexttile
    hold on
    plot(N2_offshore,-zc(1:end-1),'-','linewidth',2)
    title('N^2 offshore')
    set(gca,'fontsize',18)
    grid on
    box on
    xlabel('1/s^2')

    nexttile
    hold on
    plot(T_range_onshore,-zc,'-','linewidth',2)
    title('\Theta onshore')
    set(gca,'fontsize',18)
    grid on
    box on
    ylabel ('Depth (m)')
    xlabel('Deg C')

    nexttile
    hold on
    plot(squeeze(ambient_S_2D(1,:)),-zc,'-','linewidth',2)
    title('SP onshore')
    set(gca,'fontsize',18)
    grid on
    box on
    xlabel('PSU')

    nexttile
    hold on
    plot(dens_range_onshore,-zc,'-','linewidth',2)
    title('\rho_\theta onshore')
    set(gca,'fontsize',18)
    grid on
    box on
    xlabel('kg/m^3')

    nexttile
    hold on
    plot(N2_onshore,-zc(1:end-1),'-','linewidth',2)
    title('N^2 onshore')
    set(gca,'fontsize',18)
    grid on
    box on
    xlabel('1/s^2')

    saveas(gcf,'offshore_onshore_prof.png')
end

%% generate ptracers file and mask for CDW offshore
ptracers_file_CDW = ones(nx,ny,nz);
ptracers_file_CDW(init_T_file<1.0)=0;

% mask file is the same, same max restoration as other fields, no ramps
rbcs_mask_file_CDW = 1/24/max_restore*ptracers_file_CDW;

fid=fopen('ptracers_file_CDW.bin','w','b'); fwrite(fid,ptracers_file_CDW,acc);fclose(fid);
fid=fopen('rbcs_mask_file_CDW.bin','w','b'); fwrite(fid,rbcs_mask_file_CDW,acc);fclose(fid);

if local_plot
    % plot restoration region
    figure('pos',[200,200,700,600]);
    width=7; % relative to colorbar
    t1=tiledlayout(1,width);
    nexttile(t1,1,[1 width-1]);

    % yzdat=squeeze(rbcs_mask_file_CDW(nx/2,:,:))';
    % yzdat = 1./yzdat/24;
    % yzdat(isinf(yzdat))=100;
    % levels=[1:1:10];

    yzdat=squeeze(ptracers_file_CDW(nx/2,:,:))';
    levels=[0:0.1:1];
    numLevels=length(levels);
    cmap=flipud(cmocean('matter',numLevels));


    [latc_grid_yz,Zc_grid_yz]=meshgrid(latc,-zc);
    contourf(latc_grid_yz,Zc_grid_yz,yzdat,levels)
    hold on
    plot(latc,Ht,'-k','linewidth',2)
    plot(latc,ha,'--k','linewidth',2)
    title('CDW')
    set(gca,'fontsize',20)
    box on

    add_discrete_colorbar(t1,width,levels,cmap,'concentration',2,2)

    saveas(gcf,'cdw_region.png')
end



%% Check stability
if local_plot
    check_TS_stability3
end

%% write to plotting_data.mat
if local_plot
save('plotting_data.mat','bathy_combined','icetopo','lonc','latc','vwindkat','zc','trough_width_deg','dxfile','dyfile') 
copyfile('plotting_data.mat', '/Users/galenwilcox/Desktop/walker_zhang_research/shelfbreak1/MITgcm/shelfbreak_processing/plotting_data.mat')
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

