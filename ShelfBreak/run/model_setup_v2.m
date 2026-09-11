% plot model setup
close all

addpath('/Users/galenwilcox/Desktop/walker_zhang_research/matlab_tools/')
load plotting_data.mat

timeslice=363;
yslice=shelfbreak_ind;
xslice=nx/2;
X=lonc;
Y=latc;
[latc_grid,lonc_grid]=meshgrid(Y,X);
[lonc_grid_xz,Zc_grid_xz]=meshgrid(X,-zc);
[lonc_grid_yz,Zc_grid_yz]=meshgrid(Y,-zc);
[bathy_gridxz,Zgrid_bathyxz]=meshgrid(abs(bathy_combined(:,yslice)),zc);
[bathy_gridyz,Zgrid_bathyyz]=meshgrid(abs(bathy_combined(xslice,:)),zc);
[icetopo_gridyz,Zgrid_icetopoyz]=meshgrid(abs(icetopo(xslice,:)),zc);
zgg=zgp1(2:end);


temp_levels=[-2.0:0.25:1.5];
salt_levels=[33.8:0.05:35.0];   
%isopyc_levels=dens_top:0.1:dens_bottom;

isopyc_levels=[27.6,27.8:0.1:28.1];

%% create figure
figure;
set(gcf,'pos',[0,0,1400,500])

t = tiledlayout(1,2);
t.TileSpacing = 'compact';       % or 'compact' or 'none'
t.Padding = 'compact';

%% yz section with temperature

width=7; % relative to colorbar
t1=tiledlayout(t,1,width);
t1.Layout.Tile = 1;
nexttile(t1,1,[1 width-1]);
hold on

yzdat=squeeze(init_T_file(xslice,:,:))';
yzdat(Zgrid_bathyyz>bathy_gridyz)=NaN;
yzdat(Zgrid_icetopoyz<=icetopo_gridyz)=NaN;

levels=temp_levels;
numLevels=length(levels);
cmap=cmocean('thermal',numLevels);
yzdat(yzdat>levels(end)) = levels(end);
yzdat(yzdat<=levels(1))  = levels(1)+1e-5;


contourf(lonc_grid_yz,Zc_grid_yz,yzdat,levels,'LineColor', 'none')

% add density contours
hold on
densdat_yz=zeros(nz,ny);
for j=1:ny
    for k=1:nz
        densdat_yz(k,j)=densjmd95(squeeze(init_S_file(xslice,j,k))',squeeze(init_T_file(xslice,j,k))',0.1);
    end
end
densdat_yz=densdat_yz-1000;
densdat_yz(Zgrid_bathyyz>bathy_gridyz)=NaN;
densdat_yz(Zgrid_icetopoyz<=icetopo_gridyz)=NaN;


[C3,h3]=contour(lonc_grid_yz,Zc_grid_yz, densdat_yz,isopyc_levels,'linewidth',2,'EdgeColor','w');
clabel(C3,h3,'fontsize',14,'color','w','labelspacing',120)

ice_line=[-800,icetopo(slice,icetopo(slice,:)~=0),0];
latc_ice=latc(1:length(ice_line));

% shade area above ice shelf
yl = ylim;
x_ice = [latc_ice, fliplr(latc_ice)];
y_ice = [ice_line, repmat(yl(2), 1, length(ice_line))];
fill(x_ice, y_ice,[0.5,0.5,1],'EdgeColor', 'none', 'FaceAlpha', 0.6);
y2=plot(latc_ice,ice_line,'-','linewidth',3,'color',[0.4,0.4,1]);
text(-75.95,-70,'Ice shelf','fontsize',14,'color',[0.4,0.4,1])

% Shade area below bathy
yl = ylim;
x_bathy = [latc, fliplr(latc)];
y_bathy = [Ht, repmat(yl(1), 1, length(Ht))];
fill(x_bathy, y_bathy, [220 220 220]/255, 'EdgeColor', 'none', 'FaceAlpha', 1);
y3=plot(latc,Ht,'-k','linewidth',3);

y1=plot(latc,ha,'--k','linewidth',3);

% add a specific range limit, 20200525
xlim([latc(1) latc(end)]);
ylim([-2000,0])

title('Potential temp.')
xlabel('Latitude (\circS)')
ylabel('Depth (m)')
hold off
grid on
set(gca,'fontsize',20)
% set(gca,'YAxisLocation','right')
box on

add_discrete_colorbar(t1,width,levels,cmap,'Deg C',2,2)


%% yz section salinity

width=7; % relative to colorbar
t1=tiledlayout(t,1,width);
t1.Layout.Tile = 2;
nexttile(t1,1,[1 width-1]);
hold on

yzdat=squeeze(init_S_file(xslice,:,:))';
yzdat(Zgrid_bathyyz>bathy_gridyz)=NaN;
yzdat(Zgrid_icetopoyz<=icetopo_gridyz)=NaN;


levels=salt_levels;
numLevels=length(levels);
cmap=cmocean('haline',numLevels);

yzdat(yzdat>levels(end)) = levels(end);
yzdat(yzdat<=levels(1))  = levels(1)+1e-5;

contourf(lonc_grid_yz,Zc_grid_yz,yzdat,levels,'LineColor', 'none')

% add density contours
hold on
densdat_yz=zeros(nz,ny);
for j=1:ny
    for k=1:nz
        densdat_yz(k,j)=densjmd95(squeeze(init_S_file(xslice,j,k))',squeeze(init_T_file(xslice,j,k))',0.1);
    end
end
densdat_yz=densdat_yz-1000;
densdat_yz(Zgrid_bathyyz>bathy_gridyz)=NaN;
densdat_yz(Zgrid_icetopoyz<=icetopo_gridyz)=NaN;

[C3,h3]=contour(lonc_grid_yz,Zc_grid_yz, densdat_yz,isopyc_levels,'linewidth',2,'EdgeColor','w');
clabel(C3,h3,'fontsize',14,'color','w','labelspacing',120)

ice_line=[-800,icetopo(slice,icetopo(slice,:)~=0),0];
latc_ice=latc(1:length(ice_line));

% shade area above ice shelf
yl = ylim;
x_ice = [latc_ice, fliplr(latc_ice)];
y_ice = [ice_line, repmat(yl(2), 1, length(ice_line))];
fill(x_ice, y_ice,[0.5,0.5,1],'EdgeColor', 'none', 'FaceAlpha', 0.6);
y2=plot(latc_ice,ice_line,'-','linewidth',3,'color',[0.4,0.4,1]);
text(-75.95,-70,'Ice shelf','fontsize',14,'color',[0.4,0.4,1])

% Shade area below bathy
yl = ylim;
x_bathy = [latc, fliplr(latc)];
y_bathy = [Ht, repmat(yl(1), 1, length(Ht))];
fill(x_bathy, y_bathy, [220 220 220]/255, 'EdgeColor', 'none', 'FaceAlpha', 1);
y3=plot(latc,Ht,'-k','linewidth',3);

y1=plot(latc,ha,'--k','linewidth',3);

% add a specific range limit, 20200525
xlim([latc(1) latc(end)]);
ylim([-2000,0])

title('Practical salinity')
xlabel('Latitude (\circS)')
ylabel('Depth (m)')
hold off
grid on
set(gca,'fontsize',20)
% set(gca,'YAxisLocation','right')
box on

add_discrete_colorbar(t1,width,levels,cmap,'PSU',2,2)

exportgraphics(gcf, 'model_setup_contour_TS.png','resolution',300)

%% density

width=7; % relative to colorbar
t1=tiledlayout(t,1,width);
t1.Layout.Tile = 3;
nexttile(t1,1,[1 width-1]);
hold on

yzdat=dens_2D';
yzdat(Zgrid_bathyyz>bathy_gridyz)=NaN;
yzdat(Zgrid_icetopoyz<=icetopo_gridyz)=NaN;


levels=[27.2:0.05:28.2];
numLevels=length(levels);
cmap=cmocean('dens',numLevels);

yzdat(yzdat>levels(end)) = levels(end);
yzdat(yzdat<=levels(1))  = levels(1);

contourf(lonc_grid_yz,Zc_grid_yz,yzdat,levels,'LineColor', 'none')

% add density contours
hold on
densdat_yz=zeros(nz,ny);
for j=1:ny
    for k=1:nz
        densdat_yz(k,j)=densjmd95(squeeze(init_S_file(xslice,j,k))',squeeze(init_T_file(xslice,j,k))',0.1);
    end
end
densdat_yz=densdat_yz-1000;
densdat_yz(Zgrid_bathyyz>bathy_gridyz)=NaN;
densdat_yz(Zgrid_icetopoyz<=icetopo_gridyz)=NaN;

[C3,h3]=contour(lonc_grid_yz,Zc_grid_yz, densdat_yz,isopyc_levels,'linewidth',2,'EdgeColor','w');
clabel(C3,h3,'fontsize',14,'color','w','labelspacing',120)

ice_line=[-800,icetopo(slice,icetopo(slice,:)~=0),0];
latc_ice=latc(1:length(ice_line));

% shade area above ice shelf
yl = ylim;
x_ice = [latc_ice, fliplr(latc_ice)];
y_ice = [ice_line, repmat(yl(2), 1, length(ice_line))];
fill(x_ice, y_ice,[0.5,0.5,1],'EdgeColor', 'none', 'FaceAlpha', 0.6);
y2=plot(latc_ice,ice_line,'-','linewidth',3,'color',[0.4,0.4,1]);
text(-75.95,-70,'Ice shelf','fontsize',14,'color',[0.4,0.4,1])

% Shade area below bathy
yl = ylim;
x_bathy = [latc, fliplr(latc)];
y_bathy = [Ht, repmat(yl(1), 1, length(Ht))];
fill(x_bathy, y_bathy, [220 220 220]/255, 'EdgeColor', 'none', 'FaceAlpha', 1);
y3=plot(latc,Ht,'-k','linewidth',3);

y1=plot(latc,ha,'--k','linewidth',3);

% add a specific range limit, 20200525
xlim([latc(1) latc(end)]);
ylim([-2000,0])

title('Potential density')
xlabel('Latitude (\circS)')
ylabel('Depth (m)')
hold off
grid on
set(gca,'fontsize',20)
% set(gca,'YAxisLocation','right')
box on

add_discrete_colorbar(t1,width,levels,cmap,'kg/m^3',1,2)


%% thermal wind velocity
% integrate from z=0
% Constants
g = 9.81;          % m/s^2
rho0 = rhoConst;       % reference density (kg/m^3)
f = 2*7.29e-5*sind(-76);          % Coriolis parameter (1/s)

[~,drho_dy] = gradient(dens_2D);
drho_dy=drho_dy/(dyfile(1)*111e3);

Uvel = zeros(ny,nz);
for k=2:nz
    Uvel(:,k) = Uvel(:,k-1) + -g/(f*rho0)*drho_dy(:,k)*dzfile(k);
end


% plot it
width=7; % relative to colorbar
t1=tiledlayout(t,1,width);
t1.Layout.Tile = 4;
nexttile(t1,1,[1 width-1]);
hold on

yzdat=Uvel';
yzdat(Zgrid_bathyyz>bathy_gridyz)=NaN;
yzdat(Zgrid_icetopoyz<=icetopo_gridyz)=NaN;


levels=[-0.25:0.025:0.25];
numLevels=length(levels);
cmap=cmocean('balance',numLevels);

yzdat(yzdat>levels(end)) = levels(end);
yzdat(yzdat<=levels(1))  = levels(1);

contourf(lonc_grid_yz,Zc_grid_yz,yzdat,levels,'LineColor', 'none')

% add density contours
hold on
densdat_yz=zeros(nz,ny);
for j=1:ny
    for k=1:nz
        densdat_yz(k,j)=densjmd95(squeeze(init_S_file(xslice,j,k))',squeeze(init_T_file(xslice,j,k))',0.1);
    end
end
densdat_yz=densdat_yz-1000;
densdat_yz(Zgrid_bathyyz>bathy_gridyz)=NaN;
densdat_yz(Zgrid_icetopoyz<=icetopo_gridyz)=NaN;

[C3,h3]=contour(lonc_grid_yz,Zc_grid_yz, densdat_yz,isopyc_levels,'linewidth',2,'EdgeColor','w');
clabel(C3,h3,'fontsize',14,'color','w','labelspacing',120)

ice_line=[-800,icetopo(slice,icetopo(slice,:)~=0),0];
latc_ice=latc(1:length(ice_line));

% shade area above ice shelf
yl = ylim;
x_ice = [latc_ice, fliplr(latc_ice)];
y_ice = [ice_line, repmat(yl(2), 1, length(ice_line))];
fill(x_ice, y_ice,[0.5,0.5,1],'EdgeColor', 'none', 'FaceAlpha', 0.6);
y2=plot(latc_ice,ice_line,'-','linewidth',3,'color',[0.4,0.4,1]);
text(-75.95,-70,'Ice shelf','fontsize',14,'color',[0.4,0.4,1])

% Shade area below bathy
yl = ylim;
x_bathy = [latc, fliplr(latc)];
y_bathy = [Ht, repmat(yl(1), 1, length(Ht))];
fill(x_bathy, y_bathy, [220 220 220]/255, 'EdgeColor', 'none', 'FaceAlpha', 1);
y3=plot(latc,Ht,'-k','linewidth',3);

y1=plot(latc,ha,'--k','linewidth',3);

% add a specific range limit, 20200525
xlim([latc(1) latc(end)]);
ylim([-2000,0])

title('Geostrophic U velocity')
xlabel('Latitude (\circS)')
ylabel('Depth (m)')
hold off
grid on
set(gca,'fontsize',20)
% set(gca,'YAxisLocation','right')
box on


add_discrete_colorbar(t1,width,levels,cmap,'m/s',2,2)

exportgraphics(gcf, 'model_setup_contour.png')

%% second figure for wind and air temp
figure;
set(gcf,'pos',[800,500,1000,600])
tiledlayout(2,1)

nexttile
hold on
plot(time_6hr/4-1,vwind_kat_max,'-','linewidth',1.5)
% plot(time_6hr/4-1,vwind_rest,'-','linewidth',1.5)
plot(time_6hr/4-1,uwind_rest,'-','linewidth',1.5)
xlabel('Time (days)','fontsize',14);
axis tight
ylabel('Wind speed (m/s)','fontsize',14);
title('Wind forcing')
grid on
box on
legend('Vkat','U')
set(gca,'fontsize',20)

nexttile
yyaxis left
hold on
plot(0:15:360,airtemp_val,'-','linewidth',1.5)
ylim([-30,0])
ylabel('Deg C','fontsize',14);

yyaxis right
hold on
plot(0:15:360,lwdownfield_val,'-','linewidth',1.5)
plot(0:15:360,swdownfield_val,'-','linewidth',1.5)
xlabel('Time (days)','fontsize',14);
axis tight
ylabel('W/m^2','fontsize',14);
title('Radiative forcing')
grid on
box on
legend('Air temp','LW rad','SW rad')
set(gca,'fontsize',20)

exportgraphics(gcf, 'model_forcing_variability.png')


%% 3D bathy
figure;
set(gcf,'pos',[0,0,800,600])

width=14; % relative to colorbar
t1=tiledlayout(1,width);
nexttile(t1,1,[1 width-1]);
hold on

% Plot bathymetry
xydat=bathy_combined;
levels=[-2000:100:0];
numLevels=length(levels);
cmap=flipud(cmocean('deep',numLevels));

surf(lonc_grid, latc_grid, bathy_combined,'FaceColor', 'interp','EdgeColor', 'none');
colormap(cmap)
hold on
contour3(lonc_grid, latc_grid, bathy_combined, [-2200,0],'LineColor', 'k','LineWidth', 1);
hold on

text(lonc(round(7*nx/8)), latc(iceshelf_north_end_index-30), 70, ...
    'Coast', ...
    'Color', 'k', ...
    'FontSize', 18, ...
    'HorizontalAlignment', 'center', ...
    'Rotation', 0);


% add ice
ice_mask=icetopo~=0;
icetopo2=icetopo;
icetopo2(~ice_mask) = NaN;  % mask zeros
left_ind=find(~isnan(icetopo2(:,iceshelf_north_end_index-1)),1,'first');
right_ind=find(~isnan(icetopo2(:,iceshelf_north_end_index-1)),1,'last');
icetopo2(left_ind:right_ind,iceshelf_north_end_index)=0;

ice_levels=[-4000:200:0];
numLevels=length(ice_levels);
cmap_ice=cmocean('ice',numLevels);

ice_idx = discretize(icetopo2-0.01, ice_levels);

ice_idx_filled = ice_idx;
ice_idx_filled(isnan(ice_idx)) = 1;  % temporary index just to allow RGB mapping

% Map to RGB
[M, N] = size(icetopo2);
C_ice_rgb = reshape(cmap_ice(ice_idx_filled, :), [M, N, 3]);
C_ice_rgb(repmat(~ice_mask, [1 1 3])) = 1;  % white background

surf(lonc_grid, latc_grid, icetopo2, ...
    'CData', C_ice_rgb, ...
    'FaceColor', 'texturemap', ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 'texturemap', ...
    'AlphaData', double(ice_mask), ...
    'AlphaDataMapping', 'none');

ice_top_z = 0;

% Set non-ice regions to NaN so they're skipped
ice_top_z_data = zeros(size(icetopo2));
ice_top_z_data(isnan(icetopo2)) = NaN;

% Choose top surface color (e.g., light blue)
top_color = cmap_ice(end,:);

% Plot the flat top surface
surf(lonc_grid, latc_grid, ice_top_z_data, ...
    'FaceColor', top_color, ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 1);


hold on

text(lonc(nx/2), latc(iceshelf_north_end_index-30), 70, ...
    'Ice shelf', ...
    'Color', [0.5,0.5,1], ...
    'FontSize', 18, ...
    'HorizontalAlignment', 'center', ...
    'Rotation', 0);

% Plot wind contours on z=0 surface
vwind_to_plot = vwind(:,:,timeslice);
vwind_to_plot(bathy_combined==0) = NaN;
vwind_to_plot(:,1:iceshelf_north_end_index) = NaN;
vwind_to_plot(vwind_to_plot==0)=-0.000001;

% Compute contours manually
C = contourc(lonc, latc, vwind_to_plot', [5 10 20]);
% Parse and plot manually at z = 0
hold on
i = 1;
while i < size(C,2)
    level = C(1,i);
    npoints = C(2,i);
    x = C(1,i+1:i+npoints);
    y = C(2,i+1:i+npoints);
    z = zeros(size(x)); 
    plot3(x, y, z, '-r', 'LineWidth', 1.5);
    i = i + npoints + 1;
end
hold on
lon_wind_xaxis = lonc(vwind_to_plot(:,iceshelf_north_end_index+1) >= 1);
lat_wind = latc(iceshelf_north_end_index) * ones(size(lon_wind_xaxis));
zline = zeros(size(lon_wind_xaxis));
plot3(lon_wind_xaxis, lat_wind, zline, '-r', 'LineWidth', 1.5);

text(lonc(nx/2), latc(iceshelf_north_end_index+50), 50, ...
    'Vwind', ...
    'Color', 'r', ...
    'FontSize', 18, ...
    'HorizontalAlignment', 'center', ...
    'Rotation', 0);

% fill in bottom and add wireframe
% Choose bottom value
z_bottom = -zgp1(end);

% Extract corners
X = lonc_grid;
Y = latc_grid;
Z = bathy_combined;

% Left edge (first column)
fill3(...
    [X(:,1); flipud(X(:,1))], ...
    [Y(:,1); flipud(Y(:,1))], ...
    [Z(:,1); z_bottom*ones(size(Z(:,1)))], ...
    [0.5 0.5 0.5], 'FaceAlpha', 1, 'EdgeColor', 'k','linewidth',1);

% Right edge (last column)
fill3(...
    [X(:,end); flipud(X(:,end))], ...
    [Y(:,end); flipud(Y(:,end))], ...
    [Z(:,end); z_bottom*ones(size(Z(:,end)))], ...
    [0.5 0.5 0.5], 'FaceAlpha', 1, 'EdgeColor', 'k','linewidth',1);

% Bottom edge (first row)
fill3(...
    [X(1,:) fliplr(X(1,:))], ...
    [Y(1,:) fliplr(Y(1,:))], ...
    [Z(1,:) z_bottom*ones(size(Z(1,:)))], ...
    [0.5 0.5 0.5], 'FaceAlpha', 1, 'EdgeColor', 'k','linewidth',1);

% Top edge (last row)
fill3(...
    [X(end,:) fliplr(X(end,:))], ...
    [Y(end,:) fliplr(Y(end,:))], ...
    [Z(end,:) z_bottom*ones(size(Z(end,:)))], ...
    [0.5 0.5 0.5], 'FaceAlpha', 1, 'EdgeColor', 'k','linewidth',1);

% Define domain bounds
xlim_ = [min(lonc_grid(:)), max(lonc_grid(:))];
ylim_ = [min(latc_grid(:)), max(latc_grid(:))];
zlim_ = [-zgp1(end), 0];  % or whatever you're using

% Define the 8 corner vertices of the box
[x1, x2] = deal(xlim_(1), xlim_(2));
[y1, y2] = deal(ylim_(1), ylim_(2));
[z1, z2] = deal(zlim_(1), zlim_(2));

% Corners
corners = [
    x1 y1 z1;
    x2 y1 z1;
    x2 y2 z1;
    x1 y2 z1;
    x1 y1 z2;
    x2 y1 z2;
    x2 y2 z2;
    x1 y2 z2
];

% List of edges as pairs of corner indices
edges = [
    1 2; 2 3; 3 4; 4 1;  % bottom
    5 6; 6 7; 7 8; 8 5;  % top
    1 5; 2 6; 3 7; 4 8   % verticals
];

% Plot the edges
hold on
for i = 1:size(edges,1)
    pts = corners(edges(i,:), :);
    plot3(pts(:,1), pts(:,2), pts(:,3), 'k-', 'LineWidth', 0.8);
end

xlabel('Longitude (\circE)')
ylabel('Latitude (\circS)')
zlabel('Depth (m)')
title('Model Setup')
zlim([-zgp1(end),0])
view([135 30])
axis tight
grid on
box on
set(gca, 'FontSize', 20)

add_discrete_colorbar(t1,width,levels,cmap,'Depth (m)',0,2)

exportgraphics(gcf,'model_setup_3D.png','resolution',300)

%% plot grid and bathy
figure;
set(gcf,'pos',[200,200,800,600])

t = tiledlayout(1,1);
t.TileSpacing = 'compact';       % or 'compact' or 'none'
t.Padding = 'compact';


% bathy

width=14; % relative to colorbar
t1=tiledlayout(t,1,width);
t1.Layout.Tile = 1;
nexttile(t1,1,[1 width-1]);
hold on

xydat=bathy_combined;


levels=[-2000:100:0];
% xydat(xydat<min(levels))=min(levels);

numLevels=length(levels);
cmap=flipud(cmocean('deep',numLevels));

m_proj('lambert','lon',[X(1),X(end)],'lat',[Y(1), Y(end)]);
m_contourf(lonc_grid, latc_grid, xydat, levels, 'linewidth', 1);
title('Model grid');

% add grid lines
% --- Draw grid lines ---
% Vertical (constant i, vary j)
for i = 1:4:nx
    m_line(lonc_grid(i,:), latc_grid(i,:), 'color', [1.0 0.3 0.1], 'linewidth', 1);
end

% Horizontal (constant j, vary i)
for j = 1:4:ny
    m_line(lonc_grid(:,j), latc_grid(:,j), 'color', [1.0 0.3 0.1], 'linewidth', 1);
end



m_grid('xtick',5,'fontsize',20,'box','tickdir','out');
set(gca,'fontsize',20)
box on


add_discrete_colorbar(t1,width,levels,cmap,'Depth (m)',0,2)

exportgraphics(gcf,'grid.png','resolution',300)


%% wind forcing
figure;
set(gcf,'pos',[200,200,800,600])

timeslice=363;
width=14; % relative to colorbar
t1=tiledlayout(1,width);
nexttile(t1,1,[1 width-1]);
hold on

title('Wind forcing')
interval_x=20;
interval_y=20;
scale_factor=0.05;

hold on
%levels=[round(min(vwind,[],'all')):1:round(max(vwind,[],'all'))];
levels=[0:2:vwind_kat(1)+vwind_main(1)];
numLevels=length(levels);
cmap=flipud(cmocean('speed',numLevels));
m_proj('lambert','lon',[X(1),X(end)],'lat',[Y(1), Y(end)]);
m_contourf(lonc_grid, latc_grid, vwind(:,:,timeslice),levels, 'linewidth', 1);

m_quiver(lonc_grid(1:interval_x:end,1:interval_y:end), latc_grid(1:interval_x:end,1:interval_y:end),...
    scale_factor*uwind(1:interval_x:end,1:interval_y:end,timeslice),scale_factor*vwind(1:interval_x:end,1:interval_y:end,timeslice),...
    'autoscale','off','linewidth',1.5,'color','r','MaxHeadSize', 5)


add_coast_func(lonc,latc,bathy_combined)

m_quiver(1.25, -75.5,scale_factor*10,0,'autoscale','off','color', 'b','linewidth',2,'MaxHeadSize', 5);
m_text(1, -75.7,'10 m/s','rotation',3, 'FontSize', 20,'color','b')


m_grid('xtick',5,'fontsize',20,'box','tickdir','out');


set(gca,'fontsize',20)
box on

add_discrete_colorbar(t1,width,levels,cmap,'Vwind (m/s)',0,2)


exportgraphics(gcf,'wind_field.png','resolution',300)


