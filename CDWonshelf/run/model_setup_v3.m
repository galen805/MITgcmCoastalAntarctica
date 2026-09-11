% plot model setup
close all

xslice=nx/2;
yslice=1;
timeslice=363;
X=lonc;
Y=latc;
[latc_grid,lonc_grid]=meshgrid(Y,X);
[lonc_grid_xz,Zc_grid_xz]=meshgrid(X,-zc);
[lonc_grid_yz,Zc_grid_yz]=meshgrid(Y,-zc);
[bathy_gridxz,Zgrid_bathyxz]=meshgrid(abs(bathy_combined(:,yslice)),zc);
[bathy_gridyz,Zgrid_bathyyz]=meshgrid(abs(bathy_combined(xslice,:)),zc);
[icetopo_gridyz,Zgrid_icetopoyz]=meshgrid(abs(icetopo(xslice,:)),zc);
zgg=zgp1(2:end);


%% one model setup figure
figure;
set(gcf,'position',[200 200 1200 800])

t = tiledlayout(2,2);
t.TileSpacing = 'compact';       % or 'compact' or 'none'
t.Padding = 'compact';

%% bathymetry first
ax1 = nexttile(t, 1);   % row 1, cols 1–2
hold on

% Plot bathymetry
xydat=bathy_combined;
levels=[-800:50:0];

surf(lonc_grid, latc_grid, bathy_combined,'FaceColor', 'interp','EdgeColor', 'none');
contour3(lonc_grid, latc_grid, bathy_combined, [-2200,0],'LineColor', 'k','LineWidth', 1);

% add ice
ice_mask=icetopo~=0;
icetopo(~ice_mask) = NaN;  % mask zeros
left_ind=find(~isnan(icetopo(:,iceshelf_north_end_index-1)),1,'first');
right_ind=find(~isnan(icetopo(:,iceshelf_north_end_index-1)),1,'last');
icetopo(left_ind:right_ind,iceshelf_north_end_index)=0;

ice_levels=[-2000:200:0];
numLevels=length(ice_levels);
cmap_ice=cmocean('ice',numLevels);

ice_idx = discretize(icetopo-0.01, ice_levels);

ice_idx_filled = ice_idx;
ice_idx_filled(isnan(ice_idx)) = 1;  % temporary index just to allow RGB mapping

% Map to RGB
[M, N] = size(icetopo);
C_ice_rgb = reshape(cmap_ice(ice_idx_filled, :), [M, N, 3]);
C_ice_rgb(repmat(~ice_mask, [1 1 3])) = 1;  % white background

surf(lonc_grid, latc_grid, icetopo, ...
    'CData', C_ice_rgb, ...
    'FaceColor', 'texturemap', ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 'texturemap', ...
    'AlphaData', double(ice_mask), ...
    'AlphaDataMapping', 'none');

ice_top_z = 0;

% Set non-ice regions to NaN so they're skipped
ice_top_z_data = zeros(size(icetopo));
ice_top_z_data(isnan(icetopo)) = NaN;

% Choose top surface color (e.g., light blue)
top_color = cmap_ice(end,:);

% Plot the flat top surface
surf(lonc_grid, latc_grid, ice_top_z_data, ...
    'FaceColor', top_color, ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 1);


hold on

% Plot wind contours on z=0 surface
vwind_to_plot = vwind(:,:,timeslice);
vwind_to_plot(bathy_combined==0) = NaN;
vwind_to_plot(:,1:iceshelf_north_end_index) = NaN;
vwind_to_plot(vwind_to_plot==0)=-0.000001;

% Compute contours manually
C = contourc(lonc, latc, vwind_to_plot',4);
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

xlabel('\circE')
ylabel('\circS')
zlabel('Depth (m)')
zlim([-zgp1(end),0])
view([215 40])
axis tight
grid on
box on

% labeling
text(lonc(round(7*nx/8)), latc(iceshelf_north_end_index-20), 30, ...
    'Coast', ...
    'Color', 'k', ...
    'FontSize', 14, ...
    'HorizontalAlignment', 'center', ...
    'Rotation', 0);

text(lonc(nx/2-5), latc(iceshelf_north_end_index-20), 30, ...
    'Ice shelf', ...
    'Color', [0.5,0.5,1], ...
    'FontSize', 14, ...
    'HorizontalAlignment', 'center', ...
    'Rotation', 0);

text(lonc(nx/2+5), latc(iceshelf_north_end_index+70), 30, ...
    'Vwind', ...
    'Color', 'r', ...
    'FontSize', 14, ...
    'HorizontalAlignment', 'center', ...
    'Rotation', 0);

set(gca, 'FontSize',14)
title('(a) Model geometry','fontsize',14)

discreteColorbar(gca,levels, flipud(cmocean('deep')),'eastoutside','Depth (m)',2);


%% atmospheric forcing
t2 = tiledlayout(t,2,1);
t2.Layout.Tile = 2;

nexttile(t2)
hold on
plot((time_6hr/4-1),vwind_kat_max,'-','linewidth',1)
% plot((time_6hr/4-1)/360,vwind_rest,'-','linewidth',1.5)
plot((time_6hr/4-1),uwind_rest,'-','linewidth',1)
% xlabel('Day of year');
axis tight
ylabel('m/s');
grid on
box on
legend('V','U')
set(gca,'fontsize',14)
title('(b) Wind forcing','fontsize',14)

nexttile(t2)
yyaxis right
hold on
plot((0:15:numATEntries*15),lwdownfield_val,'--','linewidth',1.5)
plot((0:15:numATEntries*15),swdownfield_val,'-','linewidth',1.5)
xlabel('Day of year')
axis tight
ylabel('W/m^2');
grid on
box on

yyaxis left
hold on
plot((0:15:numATEntries*15),airtemp_val,'-','linewidth',1.5)
ylim([-40,10])
ylabel('\circC');


legend('AT','LW','SW','location','se')
set(gca,'fontsize',14)
title('(c) Radiative forcing','fontsize',14)





%% density and model cross section
load plotting_data.mat
ax3 = nexttile(t, 3);
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


levels=[26.9:0.1:28.0];
contourf(lonc_grid_yz,Zc_grid_yz,densdat_yz,levels,'LineColor', 'none')

ice_line=[-800,icetopo(xslice,icetopo(xslice,:)~=0),0];
latc_ice=latc(1:length(ice_line));

% Shade area below bathy
yl = ylim;
x_bathy = [latc, fliplr(latc)];
y_bathy = [Ht, repmat(yl(1), 1, length(Ht))];
fill(x_bathy, y_bathy, [220 220 220]/255, 'EdgeColor', 'none', 'FaceAlpha', 1);
y3=plot(latc,Ht,'-k','linewidth',3);

y1=plot(latc,ha,'--k','linewidth',3);

% shade area above ice shelf
yl = ylim;
x_ice = [latc_ice, fliplr(latc_ice)];
y_ice = [ice_line, repmat(yl(2), 1, length(ice_line))];
fill(x_ice, y_ice,[0.8 0.8 1],'EdgeColor', 'none', 'FaceAlpha', 0.6);
y2=plot(latc_ice,ice_line,'-','linewidth',3,'color',[0.4,0.4,1]);
text(-75.95,-70,'Ice','fontsize',14,'color',[0.4,0.4,1])
text(-75.95,-150,'Shelf','fontsize',14,'color',[0.4,0.4,1])

legend([y1,y3],'Ambient','Trough')

% add a specific range limit, 20200525
xlim([latc(1) latc(end)]);
ylim([-800,0])


title('(d) Potential density')
xlabel('Latitude (\circS)')
ylabel('Depth (m)')
hold off
grid on
set(gca,'fontsize',14)
% set(gca,'YAxisLocation','right')
box on

discreteColorbar(gca,levels, cmocean('dense'),'eastoutside','kg/m^3',2);
set(gca,'fontsize',14)


%% OBCS
t3 = tiledlayout(t,1,2);
t3.Layout.Tile = 4;

nexttile(t3,1)
hold on
plot(squeeze(init_T_file(xslice,1,:)),-zc,'-k','linewidth',2)
plot(squeeze(OBNtFile(xslice,:)),-zc,'-r','linewidth',2)
title('(e) Temperature')
legend('Initial','Forcing')
set(gca,'fontsize',14)
grid on
box on
axis tight
xlim([-2,2])
ylabel ('Depth (m)')
xlabel('\circC')

nexttile(t3,2)
hold on
plot(squeeze(init_S_file(xslice,1,:)),-zc,'-k','linewidth',2)
plot(squeeze(OBNsFile(xslice,:)),-zc,'-r','linewidth',2)
legend('Initial','Forcing')
title('(f) Practical salinity')
set(gca,'fontsize',14)
grid on
box on
axis tight
xlabel('PSU')


% exportgraphics(gcf, 'full_model_setup.png','resolution',300)

exportgraphics(gcf, 'model_setup.png','Resolution',500);



