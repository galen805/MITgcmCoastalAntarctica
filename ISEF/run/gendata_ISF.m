%clear
close all


% This is a matlab script that generates the input data

% model for WAP-style polynyas, warm shelf
% based on CDWmelt model
% created by gwilcox beginning 2/28/24

% revised 12/27/2024

%% initialize grid and time discretization
% CENTERED Dimensions of grid, high resolution region
nx1=200;
ny=250;
nz=40;

delta_z = 20;
dlat = 0.01;
dlon = 0.04;

acc = 'real*8';
ieee = 'b';

% time step
dt=30;

%% Grid file
dxfile = dlon*ones(1,nx1);
dyfile = dlat*ones(1,ny);
dzfile = delta_z*ones(1,nz);

% The above are the ones with uniform high resolution

% Now we add the one with nonuniform linearly growing grid
% On both west and east sides of x, and north of y

nx_west = 50;
nx_east = 50;
nx = nx1+nx_west+nx_east; % 200 + 25 + 25

dlon_max = 1.0;

west_add=dlon+[1:1:nx_west].^2*(dlon_max-dlon)/nx_west^2;
west_add=flip(west_add);
east_add=dlon+[1:1:nx_east].^2*(dlon_max-dlon)/nx_east^2; 

%west_add = linspace(dlon_max,dlon,nx_west);
%east_add = linspace(dlon,dlon_max,nx_east);

dxfile = [west_add dxfile east_add];

%% Add the west and east lengthen layers, keep the same north dimension
% symmetric, making the total nx to be 250
% 200 + 25*2

% print out the bin files
fid=fopen('dxfile.bin','w','b'); fwrite(fid,dxfile,acc);fclose(fid);
fid=fopen('dyfile.bin','w','b'); fwrite(fid,dyfile,acc);fclose(fid);
fid=fopen('dzfile.bin','w','b'); fwrite(fid,dzfile,acc);fclose(fid);


%% longtitude and latitude vectors at face g and center c
northen_end_lat = -74; % starting at -74 lat

latg = northen_end_lat+cumsum(dyfile);
latc = latg+dyfile/2;

long = cumsum(dxfile); % changed from -80 to -77
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

iceshelf_north_end_index = 100; % important, indicate land end
%
bathy_combined(1:(west_land_index), 1:iceshelf_north_end_index) = 0;
bathy_combined((east_land_index):nx, 1:iceshelf_north_end_index) = 0;
%

%% Round up the land corners

% Point land northwest
land_radius = 20;
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



%% WCDW1 continental shelf gwilcox 2/28/24


% main geometry points
% flip sign of max term to control prograde vs retrograde
% for now, prograde (-max)
Hc=500; % depth at back of cavity, m
hf=500; % depth at N edge of domain, m

% linearly sloped retrograde shelf
y_bathy_ambient=latc-latc(1);
ha=Hc+(hf-Hc)/y_bathy_ambient(end)*y_bathy_ambient;
ha=-ha;

%plot(y_bathy_ambient,ha) % will plot later

% modify bathy_combined
bathy_sloped=zeros(size(bathy_combined));
for i=1:ny
    bathy_sloped(:,i)=ha(i);
end

% add to bathy_combined
bathy_sloped(bathy_combined==0)=0;
bathy_combined=bathy_sloped;


%% add trough to ambient bathymetry 2/8/24
% modified to further generalize trough shape
% kind of following Pierre St-Laurent 2013 where trough was smoothed up to ice

flag_add_trough=1;

if flag_add_trough

    % trough shape parameters
    Hc_t=700; % depth at south end of domain
    hf_t=600; % trough depth at north end of domain
    
    % calculate center depth of trough
    % retrograde
    Ht=Hc_t+(hf_t-Hc_t)/y_bathy_ambient(end)*y_bathy_ambient;

    Ht=-Ht;

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
        side_scale=0.2; % deg lat
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

        % at shelfbreak
        shape_index=ny;
        lonc_centered=lonc-trough_center(shape_index);
        xcoord_km=lonc_centered*111.32*cosd(latc(shape_index));
        plot(xcoord_km,bathy_combined(:,shape_index),'-','linewidth',2)
        yline(ha(shape_index),'--')

        xlabel('km')
        ylabel('m')
        legend('Edge of ice shelf trough','Edge of ice shelf ambient','N end trough','N end ambient','location','southeast')
        xlim([-150,150])
        %ylim([min(trough_shape),max(trough_shape)])
        grid on
        set(gca,'fontsize',12)
        title('Trough shape')
        saveas(gcf,'trough_section.png')
    end
    
end



%% Ice shelf

% A more smooth ice line with Shapiro filter

% new way CDW27 1/9/24
H1=-700;
tp12=10; % turning point 1
H2=-500;
tp23=iceshelf_north_end_index-9; % turning point 2
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

plot(latc,ice_line)

% eliminate anywhere icetopo is deeper than bathy
% 1/23/24
bathy_combined(abs(icetopo)>abs(bathy_combined))=0;
icetopo(bathy_combined==0) = 0;

% make sure icetopo and bathy agree
shouldBe1=min(bathy_combined<=icetopo,[],'all') % should be 1 (always deeper than icetopo)
shouldBe0=max(bathy_combined==0 & icetopo~=0,[],'all') % should be 0 (no instances where bathy is zero (land) and ice is nonzero)


%% write bathy and plot

fid=fopen('bathy.box','w','b'); fwrite(fid,bathy_combined,acc);fclose(fid);
    
% check bathy
if local_plot
    
    % requires m_map
    figure;
    [lat_grid,lon_grid]=meshgrid(latc,lonc);
    hold on
    cmocean('haline')
    %m_proj('Transverse mercator','lon',[min(lonc) max(lonc)],'lat',[min(latc) max(latc)]);
    m_proj('Transverse mercator','lon',[lonc(end)/2-12, lonc(end)/2+12],'lat',[min(latc) max(latc)]);
    [C,h]=m_contourf(lon_grid, lat_grid, bathy_combined,10);
    %h=pcolor(bathy_combined');
    %set(h, 'EdgeColor', 'none');
    title('Model bathymetry')
    cb=colorbar;
    cb.Label.String = 'Depth (m)';
    cb.Label.FontSize = 12;
    cb.Label.Rotation = -90;
    m_grid('box','tickdir','out');
    set(gca,'fontsize',12)
    saveas(gcf,'bathy_contour.png')



    % % for use in poseidon
    % figure;
    % contourf(bathy_combined',20);
    % %h=pcolor(bathy_combined');
    % %set(h, 'EdgeColor', 'none');
    % colorbar;
    % %saveas(gcf,'bathy_contour.png')


    % % plot bathy surface in 3D
    % figure;
    % surf(lonc,latc,bathy_combined','facelighting','gouraud','EdgeAlpha',0.25)
    % colorbar;
end


%% write icetopo and plot

fid=fopen('icetopo.exp1','w','b'); fwrite(fid,icetopo,acc);fclose(fid);
% fid=fopen('pload.exp1','w','b'); fwrite(fid,-icetopo,acc);fclose(fid);

if local_plot

    % Plot iceshelf shape with bathy

    figure;
    hold on
    plot(latc(1:iceshelf_north_end_index),ice_line(1:iceshelf_north_end_index),'-','linewidth',4)
    plot(latc,Ht,'-k','linewidth',4)
    plot(latc,ha,'--k','linewidth',4)
    legend('ice shelf','trough','nominal bathy')
    xlabel('deg lat')
    ylabel('m')
    hold off
    grid on
    axis tight
    set(gca,'fontsize',12)
    saveas(gcf,'bathy_ice_yplot.png')


    % contour plot
    figure;
    contourf(icetopo',20)
    title('icetopo')
    %h=pcolor(icetopo');
    %set(h, 'EdgeColor', 'none');
    colorbar;
    saveas(gcf,'icetopo_xy.png')

    % surface plot
    figure;
    hold on
    surf(lonc,latc,bathy_combined'/1000,'facelighting','gouraud','EdgeAlpha',0.25)
    surf(lonc,latc(1:iceshelf_north_end_index),icetopo(:,1:iceshelf_north_end_index)'/1000,'facelighting','gouraud','EdgeAlpha',0.25)
    colorbar;

end



%% wind and air external forcings

% not used: vwind1 = 20*ones(size(bathy_combined)); % Wind all over

% Wind within rectangle, decaying linearly outwards
% vwind1 = zeros(size(bathy_combined));

%wind_main = 20
vwind_main=wind_main
wind_rest = 0;

% match x field to the ice shelf width
xwind1 = west_land_index; xwind2 = east_land_index;
% 100 grid points in y for the main field
ywind1 = (iceshelf_north_end_index); ywind2 = ywind1+50;

% vwind1(xwind1:xwind2, ywind1:ywind2) = wind_main;

wind_decay_length = 25;
wind_decay_frac = 1-(1:wind_decay_length)/wind_decay_length;

% for k = 1:length(wind_decay_frac)
%     vwind1(xwind1-k:xwind2+k, ywind2+k) = wind_decay_frac(k)*(wind_main-wind_rest)+wind_rest;
%     vwind1(xwind1-k, ywind1:ywind2+k) = wind_decay_frac(k)*(wind_main-wind_rest)+wind_rest;
%     vwind1(xwind2+k, ywind1:ywind2+k) = wind_decay_frac(k)*(wind_main-wind_rest)+wind_rest;
% end
% 
% 
% % relaxation points 3 cells each boundary for winds
% % not really used here
% vwind1(1:3,:) = 0;
% vwind1((end-2):end,:) = 0;
% vwind1(:,(end-2):end) = 0;

%% 20200322 overwrite- change the V-wind to a half-oval shape
mid_xwind12 = (xwind1+xwind2)/2;
vwind1 = zeros(size(bathy_combined));

wind_decay_frac_oval_flip = flip((1:wind_decay_length)/wind_decay_length);
for ki = length(wind_decay_frac_oval_flip):(-1):1
    for xi = 1:size(vwind1,1)
        for yi = ywind1:size(vwind1,2)
            
            if (xi-mid_xwind12)^2/(mid_xwind12-xwind1+ki)^2 + (yi-ywind1)^2/(ywind2-ywind1+ki)^2 <= 1
                vwind1(xi,yi) = wind_main*wind_decay_frac_oval_flip(ki);
            end
        end
    end
end


%% Switch between fluctuating wind and constant wind
wind_is_fluc = 1;
wind_is_const = (~wind_is_fluc);


% 41-1=40 records for wind field
% the 1st one is 0 all over; there is a linear interpolation in time
% 5 day for each period, 200 days for input

% fluctuating wind

% 12-hour period, 200 days


%% 20200316- modified with sporadic Vwind and Uwind

% if wind_is_fluc == 1
%
%     for j = 1:801
%         if ismember(j,(1:2:801))
%             vwind(:,:,j) = zeros(size(vwind1));
%         else
%             vwind(:,:,j) = vwind1*2; % fluctuating to 40 m/s
%         end
%     end
%
% end

load wind_scale_6hr_200days_sporadic_data.mat

% 20220301 - add a linear decay factor in spring time for 100 days
% the linear decay is from vwind1*2 (40 m/s max) to vwind1*0.5 (10 m/s max)

wind_time_decay_factor = [linspace(2, 2, length(1:600)),...
    linspace(2, 0.5, length(601:800)),...
    linspace(0.5, 0.5, length(801:1601))];


if wind_is_fluc == 1
    
    for j = 1:1601
        if j == 1
            vwind(:,:,j) = zeros(size(vwind1)); % first period is time
        end

        % no springtime wind (winter only) for WCDW1
        if (j >= 2) % winter time wind, max to 5 m/s
            vwind(:,:,j) = vwind1*wind_value_scale_6hr(j);
        end


        % if (j >= 2) && (j <= 601) % winter time wind, max to 40 m/s
        %     vwind(:,:,j) = vwind1*2 *wind_value_scale_6hr(j) /40*40; % fluctuating to 40 m/s
        % end
        % 
        % if (j > 601) % after Day 150, spring time wind, max to 10 m/s
        %     % spring time decay, changed on 2022.03.01. decay from 40 m/s to 10 m/s over 100 days
        %     vwind(:,:,j) = vwind1 *wind_time_decay_factor(j) *wind_value_scale_6hr(j) /40*40; % fluctuating to 10 m/s
        % end
    end
    
end

% constant wind
if wind_is_const == 1
    
    for j = 1:1601
        
        if j == 1
            vwind(:,:,j) = zeros(size(vwind1)); % first period is time ramp
        else
            vwind(:,:,j) = vwind1; % 20m/s constant wind
        end
    end
end


fid=fopen('vwind.bin','w','b'); fwrite(fid,vwind,acc);fclose(fid);


%% Uwind
add_Uwind = 1;

if add_Uwind == 1
    %uwind_value=[-10,-4]; % winter, spring (old)
    %uwind_value=[-5,-5];
    uwind_main=uwind_value(1)
    %uwind_value=[0,0];
    % fluctuating part not yet set up for seasonal wind forcings, commented out
    % wind drops impulsively at j <= 601 (day 150)
    
    % 2020/06/25- add a switch to uniform Uwind in the coastal region
    
    Uwind_coastal_uniform_flag = 1; % flag to switch

    % assemble seasonal wind
    if Uwind_coastal_uniform_flag == 1
        for i=1:2
            uwind_1d(i,:) = [linspace(0,0,100) linspace(uwind_value(i),uwind_value(i),100) linspace(uwind_value(i),uwind_value(i),250) linspace(uwind_value(i),0,30)];
        end
    end
    
    if Uwind_coastal_uniform_flag ~= 1
        for i=1:2
            % add on 01/30/2020- uwind 0to10m/s from 76S to 75S
            % uwind 10to0m/s from 72S to the northern end, 71S
            uwind_1d(i,:) = [linspace(0,0,100) linspace(0,uwind_value(i),100) linspace(uwind_value(i),uwind_value(i),250) linspace(uwind_value(i),0,30)];

        end
    end
    
    % the following code could be vectorized... but this works for now

    uwind1 = ones(size(bathy_combined,1),1)*uwind_1d(1,:); % winter
    uwind2 = ones(size(bathy_combined,1),1)*uwind_1d(2,:); % spring
    uwind1((end-2):end,:) = 0; % eastern BC is 0 for uwind, 3 grids
    uwind2((end-2):end,:) = 0; % eastern BC is 0 for uwind, 3 grids
    
    % 20210628 - add ramp of Uwind near all the north/east/west boundaries
    add_uwind_BC_ramp_flag = 1;
    Uwind_ramp_length=10;
    if add_uwind_BC_ramp_flag == 1
        uwind1 = function_boundary_ramp_inwards_2d_setup(bathy_combined,0,uwind_value(1),0,Uwind_ramp_length);
        uwind1(Uwind_ramp_length:(end-Uwind_ramp_length+1),1:iceshelf_north_end_index) = uwind_value(1);
        uwind2 = function_boundary_ramp_inwards_2d_setup(bathy_combined,0,uwind_value(2),0,Uwind_ramp_length);
        uwind2(Uwind_ramp_length:(end-Uwind_ramp_length+1),1:iceshelf_north_end_index) = uwind_value(2);
    end
    
    
    % 20200611- Uwind const or fluc switch, both 6-hour interval
    
    % if Uwind_is_fluc == 0, choose constant Uwind
    % remember to confirm the period length in data.exf
    % and also confirm data.obcs setup with uice and vice BC
    Uwind_is_fluc = 1;
    
    if Uwind_is_fluc==0
        % default: const Uwind
        for j = 1:1601
            if j == 1
                uwind(:,:,j) = uwind1*0; % initially zero
            elseif j > 1 && j<= 601 % winter, up to day 150
                uwind(:,:,j) = uwind1*1; %  uwind_value(1) m/s
            else % springtime, after day 150
                uwind(:,:,j) = uwind2*1; % uwind_value(2) m/s
            end
        end
    end
    
    if Uwind_is_fluc==1
        % fluctuating
        for j = 1:1601
            if j == 1
                uwind(:,:,j) = uwind1*0; % initially zero
            elseif j > 1 && j<= 601 % winter, up to day 150
                uwind(:,:,j) = uwind1*1*wind_value_scale_6hr(j); % fluctuating to uwind_value(1) m/s
            else % springtime, after day 150
                uwind(:,:,j) = uwind2*1*wind_value_scale_6hr(j); % fluctuating to uwind_value(2) m/s
            end
        end
    end
    
    fid=fopen('uwind.bin','w','b'); fwrite(fid,uwind,acc);fclose(fid);
    
end


%% Airtemp- varies linearly inwards in west/east/north relaxation region

airtemp_rest = -15;
airtemp_main = -15
airtemp_const_rest_length = 0;
airtemp_ramp_to_center_length = 0;
airtemp1 = function_boundary_ramp_inwards_2d_setup(bathy_combined,airtemp_rest,airtemp_main,airtemp_const_rest_length,airtemp_ramp_to_center_length);

airtemp2 = airtemp1;

airtemp(:,:,1) = airtemp1; % Day 0
airtemp(:,:,2) = airtemp1; % Day 50
airtemp(:,:,3) = airtemp1; % Day 100
airtemp(:,:,4) = airtemp1; % Day 150
% % every 50 day interval, after Day 150, spring time, changed on 2022.08.22.
% airtemp(:,:,5) = (-20-5)/2 * ones(size(airtemp1)); % Day 200
% airtemp(:,:,6) = (-5) * ones(size(airtemp1)); % Day 250
% % after Day 250
% airtemp(:,:,7) = (-20-5)/2 * ones(size(airtemp1)); % Day 300

% winter for 400 days WCDW1
airtemp(:,:,5) = airtemp1;
airtemp(:,:,6) = airtemp1;
airtemp(:,:,7) = airtemp1;


airtemp(:,:,8) = airtemp1; % Day 350
airtemp(:,:,9) = airtemp1; % Day 400

fid=fopen('airtempfield.bin','w','b'); fwrite(fid,airtemp,acc);fclose(fid);

%% Humidity field- const everywhere
humidfield_value = 0.0005;
humidfield1 = humidfield_value*ones(size(bathy_combined));
humidfield2 = humidfield1;

humidfield(:,:,1) = humidfield1;
humidfield(:,:,2) = humidfield2;

fid=fopen('humidfield.bin','w','b'); fwrite(fid,humidfield,acc);fclose(fid);

%% Longwave radiation field- linearly varies the same way as Air temp

lwdownfield_rest = 200;
lwdownfield_main = 200;
lwdownfield_const_rest_length = 0;
lwdownfield_ramp_to_center_length = 0;
% lwdownfield1 = function_boundary_ramp_inwards_2d_setup(bathy_combined,lwdownfield_rest,lwdownfield_main,lwdownfield_const_rest_length,lwdownfield_ramp_to_center_length);
lwdownfield_warmer = 250;
% lwdownfield1 = lwdownfield_main * ones(size(bathy_combined));
lwdownfield_1d = [linspace(lwdownfield_main,lwdownfield_main,100) linspace(lwdownfield_main,lwdownfield_warmer,380)];
lwdownfield1 = ones(size(bathy_combined,1),1)*lwdownfield_1d;

lwdownfield_uniform = lwdownfield_main * ones(size(bathy_combined));


lwdownfield(:,:,1) = lwdownfield_uniform; % Day 0
lwdownfield(:,:,2) = lwdownfield_uniform; % Day 50
lwdownfield(:,:,3) = lwdownfield_uniform; % Day 100
lwdownfield(:,:,4) = lwdownfield_uniform; % Day 150

% % every 50 day interval, after Day 150, spring time, changed on 2022.08.22.
% lwdownfield(:,:,5) = (lwdownfield_main+lwdownfield_warmer)/2 * ones(size(bathy_combined)); % Day 200
% lwdownfield(:,:,6) = (lwdownfield_warmer) * ones(size(bathy_combined)); % Day 250
% % after Day 250
% lwdownfield(:,:,7) = (lwdownfield_main+lwdownfield_warmer)/2 * ones(size(bathy_combined)); % Day 300

% winter condition only WCDW1
lwdownfield(:,:,5) = lwdownfield_uniform;
lwdownfield(:,:,6) = lwdownfield_uniform;
lwdownfield(:,:,7) = lwdownfield_uniform;

lwdownfield(:,:,8) = lwdownfield_uniform; % Day 350
lwdownfield(:,:,9) = lwdownfield_uniform; % Day 400

fid=fopen('lwdownfield.bin','w','b'); fwrite(fid,lwdownfield,acc);fclose(fid);


%% Shortwave radiation field- linearly varies the same way as Air temp
swdownfield_rest = 0;
swdownfield_main = 0;
swdownfield_const_rest_length = 0;
swdownfield_ramp_to_center_length = 0;
% swdownfield1 = function_boundary_ramp_inwards_2d_setup(bathy_combined,swdownfield_rest,swdownfield_main,swdownfield_const_rest_length,swdownfield_ramp_to_center_length);
swdownfield_warmer = 400;
% swdownfield1 = swdownfield_main * ones(size(bathy_combined));
swdownfield_1d = [linspace(swdownfield_main,swdownfield_main,100) linspace(swdownfield_main,swdownfield_warmer,380)];
swdownfield1 = ones(size(bathy_combined,1),1)*swdownfield_1d;

swdownfield_uniform = swdownfield_main * ones(size(bathy_combined));


swdownfield(:,:,1) = swdownfield_uniform; % Day 0
swdownfield(:,:,2) = swdownfield_uniform; % Day 50
swdownfield(:,:,3) = swdownfield_uniform; % Day 100
swdownfield(:,:,4) = swdownfield_uniform; % Day 150


% % every 50 day interval, after Day 150, spring time, changed on 2022.08.22.
% swdownfield(:,:,5) = (swdownfield_main+swdownfield_warmer)/2 * ones(size(bathy_combined)); % Day 200
% swdownfield(:,:,6) = (swdownfield_warmer) * ones(size(bathy_combined)); % Day 250
% % after Day 250
% swdownfield(:,:,7) = (swdownfield_main+swdownfield_warmer)/2 * ones(size(bathy_combined)); % Day 300

% winter only CDW27
swdownfield(:,:,5) = swdownfield_uniform;
swdownfield(:,:,6) = swdownfield_uniform;
swdownfield(:,:,7) = swdownfield_uniform;

swdownfield(:,:,8) = swdownfield_uniform; % Day 350
swdownfield(:,:,9) = swdownfield_uniform; % Day 400

fid=fopen('swdownfield.bin','w','b'); fwrite(fid,swdownfield,acc);fclose(fid);

%% Sea ice initial files
SIarea_init_value = 1;
SIarea_init = SIarea_init_value*ones(size(bathy_combined));
SIarea_init(bathy_combined == 0) = 0; % get rid of the land
SIarea_init(icetopo ~= 0) = 0; % get rid of ice shelf, ice tongue and pack ice
fid=fopen('SIarea_init.bin','w','b'); fwrite(fid,SIarea_init,acc);fclose(fid);

SIheff_init_value = 0.5;
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


%% 2020/05/22- Switch sea ice boundary conditions
% depending on whether Uwind is on or off

% Original sea ice BC by default

%% Seaice OBCS nx/y * time levels
% ice area = 1 west east, with decay northward
% 20 points used for decaying, about 0.8deg in lat
% north point is 0
ice_area = 1;
ice_area_bc_decay_num = 10;
ice_area_bc_decay = linspace(ice_area,0,ice_area_bc_decay_num);
ice_area_bc_decay = ice_area_bc_decay';

iceOBCSawest1 = ice_area*ones(size(bathy_combined,2),1);
iceOBCSawest1((end-ice_area_bc_decay_num+1):end)=ice_area_bc_decay;
iceOBCSawest1(1:iceshelf_north_end_index,:) = 0; % eliminate the land part
iceOBCSawest2 = iceOBCSawest1;
iceOBCSawest(:,1) = iceOBCSawest1*0.0;
iceOBCSawest(:,2) = iceOBCSawest2*0.0;
fid=fopen('iceOBCSawest.bin','w','b'); fwrite(fid,iceOBCSawest,acc);fclose(fid);


iceOBCSaeast = iceOBCSawest;
fid=fopen('iceOBCSaeast.bin','w','b'); fwrite(fid,iceOBCSaeast,acc);fclose(fid);

iceOBCSanorth1 = 0.0*ones(size(bathy_combined,1),1);
iceOBCSanorth2 = iceOBCSanorth1;
iceOBCSanorth(:,1) = iceOBCSanorth1;
iceOBCSanorth(:,2) = iceOBCSanorth2;
fid=fopen('iceOBCSanorth.bin','w','b'); fwrite(fid,iceOBCSanorth,acc);fclose(fid);


% ice thickness
% enforcing the same pattern as ice area bc
ice_thickness=0.5;

iceOBCShwest = 0.0*ice_thickness*iceOBCSawest;
fid=fopen('iceOBCShwest.bin','w','b'); fwrite(fid,iceOBCShwest,acc);fclose(fid);

iceOBCSheast = 0.0*ice_thickness*iceOBCSaeast;
fid=fopen('iceOBCSheast.bin','w','b'); fwrite(fid,iceOBCSheast,acc);fclose(fid);

iceOBCShnorth = 0.0*iceOBCSanorth;
fid=fopen('iceOBCShnorth.bin','w','b'); fwrite(fid,iceOBCShnorth,acc);fclose(fid);


%% with Uwind condition, change sea ice BC
if add_Uwind == 1   
    % add on 01/22- add sea ice u-velocity on the east/west boundaries
    % velocity is 0.0m/s
    uice_velocity=0.0;
    uiceOBCSwest = uice_velocity*ones(size(iceOBCSawest));
    uiceOBCSwest(1:iceshelf_north_end_index,:) = 0;
    fid=fopen('uiceOBCSwest.bin','w','b'); fwrite(fid,uiceOBCSwest,acc);fclose(fid);
    
    uiceOBCSeast = uiceOBCSwest;
%     % add on 01/26- Version3, use previous output from 01/25 for east BC Uice
%     load uice_xin_353.mat
%     uiceOBCSeast = [uice_xin_353' uice_xin_353'];
    
    % add on 01/29 enforce east BC Uice=0
    uiceOBCSeast = 0.0*uiceOBCSeast;
    fid=fopen('uiceOBCSeast.bin','w','b'); fwrite(fid,uiceOBCSeast,acc);fclose(fid);
    
%     %% add on 01/25- change the east sea ice thickness to be 0
%     iceOBCSheast = 0.0*iceOBCSaeast;
%     fid=fopen('iceOBCSheast.bin','w','b'); fwrite(fid,iceOBCSheast,acc);fclose(fid);
    
    
    % add on 20200213- if at boundaries the pack ice with icetopo~=0,
    % then sea ice BC=0
    % rewrite and regenerate west sea ice BC
    icetopo_at_west_BC(:,1) = icetopo(1,:);
    icetopo_at_west_BC(:,2) = icetopo(1,:);
    
    iceOBCSawest(icetopo_at_west_BC~=0) = 0;
    fid=fopen('iceOBCSawest.bin','w','b'); fwrite(fid,iceOBCSawest,acc);fclose(fid);
    
    iceOBCShwest(icetopo_at_west_BC~=0) = 0;
    fid=fopen('iceOBCShwest.bin','w','b'); fwrite(fid,iceOBCShwest,acc);fclose(fid);
    
    uiceOBCSwest(icetopo_at_west_BC~=0) = 0;
    fid=fopen('uiceOBCSwest.bin','w','b'); fwrite(fid,uiceOBCSwest,acc);fclose(fid);
    
end


%% Ptracers for meltwater and salt (tied to shelfice and seaice packages)
% Add the parts for ptracers surface intial file

ptracers2d=ones(size(bathy_combined));
ptracers2d(:,1:iceshelf_north_end_index) = 0;

ptracers3d = zeros(size(bathy_combined,1),size(bathy_combined,2),nz);

% for inp = 1:1 % here initialized with 1 all over the surface
%
%     ptracers3d(:,:,inp) = ptracers2d;
%
% end

fid=fopen('ptracers_file.bin','w',ieee); fwrite(fid,ptracers3d,acc); fclose(fid);

% Add mask file for Ptracers, use the below for RBCS

% The west/east/north boundary cells relaxation for ptracers
ptracers_mask = zeros(size(ptracers3d));

ptr_rbcs_ramp_length = 10;

ptracers_boundary_value = 1/24; % 24*3600s 1-daily relaxation

ptracers2d_bc = function_boundary_ramp_inwards_2d_setup(bathy_combined,ptracers_boundary_value,0,0,ptr_rbcs_ramp_length);

for depth_in = 1:nz
    ptracers_mask(:,:,depth_in) = ptracers2d_bc;
end

% The surface relaxation for constant supply of ptracers (to mimic salt)
ptracers_mask(ptracers3d~=0) = 1;

fid=fopen('rbcs_mask_ptracers.bin','w',ieee); fwrite(fid,ptracers_mask,acc); fclose(fid);

%plot_field(ptracers_mask,length(zc)-5,245,lonc,latc,zc,bathy_combined,icetopo)
%plot_field(ptracers3d,length(zc)-5,245,lonc,latc,zc,bathy_combined,icetopo)

%% generate salinity stratification (1D file)
%salt_top=32.9;
%salt_top=33.5;
salt_top
salt_bottom=34.5

if linstrat==1
    % linear stratification down to bottom of trough
    strat_stop=find(zc>=Hc_t,1)-2; %vertical steps
    salt_range=salt_top+(salt_bottom-salt_top)/strat_stop*[0:strat_stop-1];
    salt_range=[salt_range,salt_bottom*ones(1,nz-strat_stop)];
else
    % mixed layer form
    ML_stop=5; %100m mixed layer
    salt_stop=20; % stop at 500m
    
    salt_range=salt_top+(salt_bottom-salt_top)/(salt_stop-ML_stop)*[1:salt_stop-ML_stop];
    salt_range=[salt_top*ones(1,ML_stop),salt_range,salt_bottom*ones(1,nz-salt_stop)];

end

if local_plot

    xslice=nx/2;
    yzdat=repmat(salt_range',1,length(latc));
    [bathy_grid,Zgrid_bathy]=meshgrid(abs(bathy_combined(xslice,:)),zc);
    yzdat(bathy_grid<=Zgrid_bathy)=NaN;

    figure
    hold on
    %contourf(lat_grid_side, Zgrid, yzdat,'LineColor', 'none');
    h = pcolor(latc,-zc,yzdat);
    plot(latc,bathy_combined(xslice,:),'-k','linewidth',1.5)
    plot(latc(1:iceshelf_north_end_index),icetopo(xslice,1:iceshelf_north_end_index),'--k','linewidth',1.5)
    legend('property','bathy','iceshelf')
    set(h, 'EdgeColor', 'none');
    axis tight
    grid on
    colorbar;
    title('Salt stratification')
    saveas(gcf,'S_stratification_plot.png')

end


% write bin file
fid=fopen('salt_init_file.bin','w','b'); fwrite(fid,salt_range,acc);fclose(fid);

% print text
salt_txt = sprintf('%.4f, ', salt_range)


%% generate temperature stratification (1D file)
% similar to salinity stratification
% kind of following Pierre St-Laurent 2013
T_top=-1.8;
T_bottom=1;

ML_stop=5; %100m mixed layer
T_stop=20; % stop at 400m

if linstrat==1
    % linear stratification down to bottom of trough
    strat_stop=find(zc>=Hc_t,1)-2; %vertical steps
    T_range=T_top+(T_bottom-T_top)/strat_stop*[0:strat_stop-1];
    T_range=[T_range,T_bottom*ones(1,nz-strat_stop)];
else
    % mixed layer form
    ML_stop=5; %100m mixed layer
    strat_stop=20; % stop at 400m
    
    T_range=T_top+(T_bottom-T_top)/(strat_stop-ML_stop)*[1:strat_stop-ML_stop];
    T_range=[T_top*ones(1,ML_stop),T_range,T_bottom*ones(1,nz-strat_stop)];

end

if local_plot

    xslice=nx/2;
    yzdat=repmat(T_range',1,length(latc));
    [bathy_grid,Zgrid_bathy]=meshgrid(abs(bathy_combined(xslice,:)),zc);
    yzdat(bathy_grid<=Zgrid_bathy)=NaN;

    figure
    hold on
    %contourf(lat_grid_side, Zgrid, yzdat,'LineColor', 'none');
    h = pcolor(latc,-zc,yzdat);
    plot(latc,bathy_combined(xslice,:),'-k','linewidth',1.5)
    plot(latc(1:iceshelf_north_end_index),icetopo(xslice,1:iceshelf_north_end_index),'--k','linewidth',1.5)
    legend('property','bathy','iceshelf')
    set(h, 'EdgeColor', 'none');
    axis tight
    grid on
    colorbar;
    title('Temp stratification')
    saveas(gcf,'T_stratification_plot.png')

end


mean_T_below_300m=mean(T_range(16:end))

% write bin file
fid=fopen('T_init_file.bin','w','b'); fwrite(fid,T_range,acc);fclose(fid);

% print text
T_txt = sprintf('%.4f, ', T_range)


%% 2/24/24 gwilcox CDW temp restoration in east side of trough

% relax field files
ambient_T_2D=repmat(T_range',1,length(latc))';
for i=1:length(lonc)
    ambient_temp(i,:,:)=ambient_T_2D;
end


CDW_temp=1.5; % quite warm

% define broad limits of restoration
xlowerlim=95+nx_east;
xupperlim=170+nx_east;
ylowerlim=200;
zlowerlim=20;

rbcs_T_file=ambient_temp;
rbcs_T_file(xlowerlim:xupperlim,ylowerlim:end,zlowerlim:end)=CDW_temp;

% add ramp to top and sides
ramp_length_z=3; % grid cells
for i=1:ramp_length_z
    % upper part
    current_slice=rbcs_T_file(xlowerlim:end,ylowerlim:end,zlowerlim+i-1);
    ambient_slice=ambient_temp(xlowerlim:end,ylowerlim:end,zlowerlim+i-1);
    mz=sqrt(abs(ambient_slice-current_slice))/ramp_length_z;
    rbcs_T_file(xlowerlim:end,ylowerlim:end,zlowerlim+i-1)= current_slice  -(mz*(ramp_length_z-i)).^2;
end
rbcs_T_file(rbcs_T_file<ambient_temp)=ambient_temp(rbcs_T_file<ambient_temp);

% W side
ramp_length_x=10;
for i=1:ramp_length_x
    current_slice=rbcs_T_file(xlowerlim+i-1,ylowerlim:end,zlowerlim:end);
    ambient_slice=ambient_temp(xlowerlim+i-1,ylowerlim:end,zlowerlim:end);
    mx=sqrt(abs(ambient_slice-current_slice))/ramp_length_x;
    rbcs_T_file(xlowerlim+i-1,ylowerlim:end,zlowerlim:end)= current_slice - (mx*(ramp_length_x-i)).^2;
end
% correct corners of ramp
rbcs_T_file(rbcs_T_file<ambient_temp)=ambient_temp(rbcs_T_file<ambient_temp);

% E side
for i=1:ramp_length_x
    current_slice=rbcs_T_file(xupperlim-i+1,ylowerlim:end,zlowerlim:end);
    ambient_slice=ambient_temp(xupperlim-i+1,ylowerlim:end,zlowerlim:end);
    mx=sqrt(abs(ambient_slice-current_slice))/ramp_length_x;
    rbcs_T_file(xupperlim-i+1,ylowerlim:end,zlowerlim:end)= current_slice - (mx*(ramp_length_x-i)).^2;
end
% correct corners of ramp
rbcs_T_file(rbcs_T_file<ambient_temp)=ambient_temp(rbcs_T_file<ambient_temp);



% set up ptracer field for CDW intrusion
ptracers_CDW=rbcs_T_file-ambient_temp;
ptracers_CDW=ptracers_CDW./(CDW_temp-ambient_temp); % concentration of CDW

%% calculate mean temp in cavity
Temp=rbcs_T_file;
for i=1:nx
    for j=1:ny
        for k=1:nz
            if -zc(k)<bathy_combined(i,j)
                Temp(i,j,k,:)=NaN;
            elseif -zc(k)>icetopo(i,j)
                Temp(i,j,k,:)=NaN;
            end
        end
    end
end

T_zavg=squeeze(mean(Temp,3,'omitnan'));
T_incavity=T_zavg;
T_incavity(icetopo==0)=0;
T_avg=mean(T_incavity,'all')

%% correct density by adding a salinity anomaly
ambient_S_2D=repmat(salt_range',1,length(latc))';
for i=1:length(lonc)
    ambient_S(i,:,:)=ambient_S_2D;
end

gravity=9.81;
rhoConst = 1030;

rbcs_S_file=ambient_S;
p=0.1; % use potential density for this (surface pressure)

dens_top=densjmd95(salt_top,T_top,p);
dens_bottom=densjmd95(salt_bottom,T_bottom,p);
density_anom = 0.1*(dens_bottom-dens_top);
for i=1:nx
    for j=1:ny
        for k=1:nz
            if rbcs_T_file(i,j,k)>T_range(k)
                n_iter=0;
                rescheck=1;

                match_dens = densjmd95(salt_range(k), T_range(k), p) + density_anom;

                % Initial guess
                S = rbcs_S_file(i,j,k);

                while rescheck > 1e-5
                    % Compute current density and residual
                    dens = densjmd95(S, rbcs_T_file(i,j,k), p);
                    res = dens - match_dens;

                    % Finite difference approximation of d(dens)/dS
                    dS = 1e-4;
                    dens_dS = densjmd95(S + dS, rbcs_T_file(i,j,k), p);
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

                rbcs_S_file(i,j,k) = S;


                % my version, improved version above (from chatgpt, sped up at least 100x)
                % match_dens=densjmd95(salt_range(k),T_range(k),p) + density_anom;
                % while rescheck>1e-5
                %     rbcs_S_file(i,j,k)=rbcs_S_file(i,j,k)+1e-4;
                %     check_dens=densjmd95(rbcs_S_file(i,j,k),rbcs_T_file(i,j,k),p);
                %     res=abs(check_dens-match_dens);
                %     if res>rescheck
                %         rescheck=0;
                %         fprintf('Pick smaller increment for iteration')
                %     else
                %         rescheck=res;
                %     end
                %     n_iter=n_iter+1;
                % end
            end
        end
    end
end

%% write
fid=fopen('rbcs_S_file.bin','w','b'); fwrite(fid,rbcs_S_file,acc);fclose(fid);
fid=fopen('rbcs_T_file.bin','w','b'); fwrite(fid,rbcs_T_file,acc);fclose(fid);
fid=fopen('ptracers_file_CDW.bin','w','b'); fwrite(fid,ptracers_CDW,acc);fclose(fid);

%% mask files
% open-ocean restoration
% set mask which governs restoration (want restoration every 1 day, timescale tau set in data.rbcs)
% nonzero mask value where we want relaxation/restoration
% all restoration is the same for now

max_restore=1/24/6.5; % restoration every 6.5 days
rbcs_mask_T=zeros(size(ptracers_CDW));
rbcs_mask_T(ptracers_CDW>0)=max_restore;
%rbcs_mask_T=ptracers_CDW*max_restore;
%rbcs_mask_T(1:xlowerlim-1,1:ylowerlim-1,1:zlowerlim-1)=0; % only restore in outer sponge layer

% no ramp in restore region, want it to be like a pulse, sudden inflow every 6.5 days

ptracers_mask_CDW=rbcs_mask_T; % same restoration as temp
rbcs_mask_S=rbcs_mask_T;

% assemble rbcs mask bin files
fid=fopen('rbcs_mask_S.bin','w','b'); fwrite(fid,rbcs_mask_T,acc);fclose(fid);
fid=fopen('rbcs_mask_T.bin','w','b'); fwrite(fid,rbcs_mask_T,acc);fclose(fid);
fid=fopen('rbcs_mask_ptracers_CDW.bin','w','b'); fwrite(fid,ptracers_mask_CDW,acc);fclose(fid);

%% Check stability
if local_plot
    check_TS_stability3
end

%% write to plotting_data.mat
if local_plot
save('plotting_data.mat','bathy_combined','icetopo','lonc','latc','vwind1','zc','trough_width_deg','dxfile','dyfile') 
copyfile('plotting_data.mat', '/Users/galenwilcox/Desktop/walker_zhang_research/MITgcm/pisces/ISFwarm/ISF_processing/plotting_data.mat')
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

