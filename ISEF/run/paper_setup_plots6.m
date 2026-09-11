%
close all

yslice=245;
X=lonc;
Y=latc;
[latc_grid,lonc_grid]=meshgrid(Y,X);
[lonc_grid_xz,Zc_grid_xz]=meshgrid(X,-zc);
[bathy_gridxz,Zgrid_bathyxz]=meshgrid(abs(bathy_combined(:,yslice)),zc);

%% create figure
figure;
set(gcf,'pos',[0,0,1000,800])

t = tiledlayout(2,2);
t.TileSpacing = 'compact';       % or 'compact' or 'none'
t.Padding = 'compact';

%% bathy
width=14; % relative to colorbar
t1=tiledlayout(t,1,width);
t1.Layout.Tile = 1;
nexttile(t1,1,[1 width-1]);
hold on

levels=[-700:50:0];
numLevels=length(levels);
cmap=flipud(cmocean('deep',numLevels));

m_proj('lambert','lon',[trough_center(1)-5, trough_center(1)+5],'lat',[Y(1), Y(end)]);
m_contourf(lonc_grid, latc_grid, bathy_combined, levels, 'linewidth', 1);
title('Model bathymetry');
m_grid('xtick',4,'fontsize',20,'box','tickdir','out');


% add wind restoration
hold on
vwind_to_plot = vwind(:,:,340);
vwind_to_plot(bathy_combined==0) = NaN;
vwind_to_plot(icetopo<0) = NaN;
vwind_to_plot(vwind_to_plot==0)=-0.000001;
[C1, h1]=m_contour(lonc_grid, latc_grid, vwind_to_plot,[0,15],'--w','linewidth',1.5);
clabel(C1, h1, 'FontSize', 16, 'Color', 'w');
hold on
lon_wind_xaxis = lonc(vwind_to_plot(:,100)>=0);
m_plot(lon_wind_xaxis,latc(100)*ones(size(lon_wind_xaxis)),'--w','linewidth',1.5)

%add restoration
restore_lin=[find(squeeze(rbcs_mask_T(:,:,25))'~=0,1,"first"),find(squeeze(rbcs_mask_T(:,:,25))'~=0,1,"last")];
[restore_y,restore_x]=ind2sub(size(bathy_combined'),restore_lin);
xRect = [restore_x(1), restore_x(2), restore_x(2), restore_x(1)];
yRect = [restore_y(1), restore_y(1), restore_y(2), restore_y(2)];
m_patch(lonc(xRect), sort(latc(yRect)), 'w', 'FaceAlpha', 0.0, 'EdgeColor', 'r','linewidth',2,'linestyle','--');


% add icetopo patch
hold on
mask = icetopo ~= 0;
% Generate patches for each cell where icetopo is nonzero
for i = 1:size(lonc_grid,1)-1
    for j = 1:size(lonc_grid,2)-1
        if mask(i,j)
            x_patch = [lonc_grid(i,j), lonc_grid(i+1,j), lonc_grid(i+1,j+1), lonc_grid(i,j+1)]-0.02;
            y_patch = [latc_grid(i,j), latc_grid(i+1,j), latc_grid(i+1,j+1), latc_grid(i,j+1)];
            m_patch(x_patch, y_patch,[0.5,0.5,1], 'FaceAlpha', 0.6,'edgecolor','none');
        end
    end
end

m_grid('xtick',4,'fontsize',20,'box','tickdir','out');
set(gca,'fontsize',20)
box on

m_text(21.55,-73.5,'Ice shelf','fontsize',14,'color',[0.5,0.5,1])
m_text(21.75,-72.8,'Vwind','fontsize',14,'color','w')
m_text(22.65,-71.6,'Restoration','fontsize',14,'color','r')


add_discrete_colorbar(t1,width,levels,cmap,'Depth (m)',0,2)


%% yz section
nexttile(t)
hold on

ylim([-700,0])
ice_line=[-700,icetopo(slice,icetopo(slice,:)~=0),0];
latc_ice=latc(1:length(ice_line));

% shade area above ice shelf
yl = ylim;
x_ice = [latc_ice, fliplr(latc_ice)];
y_ice = [ice_line, repmat(yl(2), 1, length(ice_line))];
fill(x_ice, y_ice,[0.5,0.5,1],'EdgeColor', 'none', 'FaceAlpha', 0.6);
y2=plot(latc_ice,ice_line,'-','linewidth',3,'color',[0.5,0.5,1]);

% Shade area below bathy
yl = ylim;
x_bathy = [latc, fliplr(latc)];
y_bathy = [Ht, repmat(yl(1), 1, length(Ht))];
fill(x_bathy, y_bathy, [220 220 220]/255, 'EdgeColor', 'none', 'FaceAlpha', 1);
y3=plot(latc,Ht,'-k','linewidth',3);

y1=plot(latc,ha,'--k','linewidth',3);

%add restoration
restore_lin=[find(squeeze(rbcs_mask_T(slice,:,:))',1,"first"),find(squeeze(rbcs_mask_T(slice,:,:))',1,"last")];
[restore_y,restore_x]=ind2sub(size(yzdat),restore_lin);
xRect = [restore_x(1), restore_x(2), restore_x(2), restore_x(1)];
yRect = [restore_y(1), restore_y(1), restore_y(2), restore_y(2)];
y4=patch(latc(xRect), sort(-zc(yRect)), 'w', 'FaceAlpha', 0, 'EdgeColor', 'r','linewidth',2,'linestyle','--');

% add a specific range limit, 20200525
add_range_limit_flag = 1;
if add_range_limit_flag
    xlim([latc(1) latc(end)]);
    ylim([-700,0])
end

legend([y2,y3,y1,y4],'Ice shelf','Trough center','Ambient','Restoration')
title('Cross-shore section')
xlabel('Latitude (\circS)')
ylabel('Depth (m)','rotation',-90)
hold off
grid on
set(gca,'fontsize',20)
set(gca,'YAxisLocation','right')
box on


%% initial temp xz
width=7; % relative to colorbar
t1=tiledlayout(t,1,width);
t1.Layout.Tile = 3;
nexttile(t1,1,[1 width-1]);
hold on

xzdat=squeeze(rbcs_T_file(:,yslice,:))';
xzdat(Zgrid_bathyxz>bathy_gridxz)=NaN;

levels=[-2:0.25:1.5];
numLevels=length(levels);
cmap=cmocean('thermal',numLevels);

contourf(lonc_grid_xz,Zc_grid_xz,xzdat,levels,'LineColor', 'none')

% add bathy
hold on
bathy_line=bathy_combined(bathy_combined(:,yslice)~=0,yslice)';
lonc_bathy=lonc(bathy_combined(:,yslice)~=0);
yl = ylim;
x_bathy = [lonc_bathy, fliplr(lonc_bathy)];
y_bathy = [bathy_line, repmat(yl(1), 1, length(bathy_line))];
fill(x_bathy, y_bathy, [220,220,220]/255, 'EdgeColor', 'none', 'FaceAlpha', 1);
plot(lonc_bathy,bathy_line,'-k','linewidth',2); hold on;

% add restoration
restore_lin=[find(squeeze(rbcs_mask_T(:,yslice,:))',1,"first"),find(squeeze(rbcs_mask_T(:,yslice,:))',1,"last")];
[restore_y,restore_x]=ind2sub(size(xzdat),restore_lin);
xRect = [restore_x(1), restore_x(2), restore_x(2), restore_x(1)];
yRect = [restore_y(1), restore_y(1), restore_y(2), restore_y(2)];
patch(lonc(xRect), sort(-zc(yRect)), 'w', 'FaceAlpha', 0.0, 'EdgeColor', 'r','linewidth',2,'linestyle','--');

% add density contours
hold on
densdat_xz=zeros(nz,nx);
for j=1:nx
    for k=1:nz
        densdat_xz(k,j)=densjmd95(squeeze(rbcs_S_file(j,yslice,k))',squeeze(rbcs_T_file(j,yslice,k))',0.1);
    end
end
densdat_xz=densdat_xz-1000;
densdat_xz(Zgrid_bathyxz>bathy_gridxz)=NaN;
isopyc_levels = round(linspace(min(densdat_xz,[],'all'),max(densdat_xz,[],'all'),7),1);
isopyc_levels = isopyc_levels(2:end-1);
[C3,h3]=contour(lonc_grid_xz,Zc_grid_xz, densdat_xz,isopyc_levels,'linewidth',2,'EdgeColor','w');
clabel(C3,h3,'fontsize',14,'color','w','labelspacing',120)

% title and stuff
title('Potential temp -71.5 S')
ylabel('Depth (m)')
xlabel('Longitude (\circE)')
xlim([19, 26])
ylim([-700,0])
box on
set(gca,'fontsize',20)

add_discrete_colorbar(t1,width,levels,cmap,'Deg C',2,2)

%% salinity xz
width=7; % relative to colorbar
t1=tiledlayout(t,1,width);
t1.Layout.Tile = 4;
nexttile(t1,1,[1 width-1]);
hold on

xzdat=squeeze(rbcs_S_file(:,yslice,:))';
xzdat(Zgrid_bathyxz>bathy_gridxz)=NaN;

levels=linspace(min(xzdat,[],'all'),max(xzdat,[],'all'),10);
numLevels=length(levels);
cmap=cmocean('haline',numLevels);

contourf(lonc_grid_xz,Zc_grid_xz,xzdat,levels,'LineColor', 'none')

% add bathy
hold on
bathy_line=bathy_combined(bathy_combined(:,yslice)~=0,yslice)';
lonc_bathy=lonc(bathy_combined(:,yslice)~=0);
yl = ylim;
x_bathy = [lonc_bathy, fliplr(lonc_bathy)];
y_bathy = [bathy_line, repmat(yl(1), 1, length(bathy_line))];
fill(x_bathy, y_bathy, [220,220,220]/255, 'EdgeColor', 'none', 'FaceAlpha', 1);
plot(lonc_bathy,bathy_line,'-k','linewidth',2); hold on;

% add restoration
restore_lin=[find(squeeze(rbcs_mask_T(:,yslice,:))',1,"first"),find(squeeze(rbcs_mask_T(:,yslice,:))',1,"last")];
[restore_y,restore_x]=ind2sub(size(xzdat),restore_lin);
xRect = [restore_x(1), restore_x(2), restore_x(2), restore_x(1)];
yRect = [restore_y(1), restore_y(1), restore_y(2), restore_y(2)];
patch(lonc(xRect), sort(-zc(yRect)), 'w', 'FaceAlpha', 0.0, 'EdgeColor', 'r','linewidth',2,'linestyle','--');

% add density contours
hold on
densdat_xz=zeros(nz,nx);
for j=1:nx
    for k=1:nz
        densdat_xz(k,j)=densjmd95(squeeze(rbcs_S_file(j,yslice,k))',squeeze(rbcs_T_file(j,yslice,k))',0.1);
    end
end
densdat_xz=densdat_xz-1000;
densdat_xz(Zgrid_bathyxz>bathy_gridxz)=NaN;
isopyc_levels = round(linspace(min(densdat_xz,[],'all'),max(densdat_xz,[],'all'),7),1);
isopyc_levels = isopyc_levels(2:end-1);
[C3,h3]=contour(lonc_grid_xz,Zc_grid_xz, densdat_xz,isopyc_levels,'linewidth',2,'EdgeColor','w');
clabel(C3,h3,'fontsize',14,'color','w','labelspacing',120)

% title and stuff
title('Practical salinity -71.5 S')
axis tight
ylabel('Depth (m)')
xlabel('Longitude (\circE)')
xlim([19, 26])
ylim([-700,0])
box on
set(gca,'fontsize',20)

add_discrete_colorbar(t1,width,levels,cmap,'PSU',2,1)

%% annotate and save
annotation('textbox', [.06 .48,.1,.1], 'String', '(a)', 'EdgeColor', 'none', 'FontSize', 30, 'FontWeight', 'bold');
annotation('textbox', [.56 .48,.1,.1], 'String', '(b)', 'EdgeColor', 'none', 'FontSize', 30, 'FontWeight', 'bold');
annotation('textbox', [.06 0,.1,.1], 'String', '(c)', 'EdgeColor', 'none', 'FontSize', 30, 'FontWeight', 'bold');
annotation('textbox', [.56 0,.1,.1], 'String', '(d)', 'EdgeColor', 'none', 'FontSize', 30, 'FontWeight', 'bold');


exportgraphics(gcf, 'Poster_model_setup.png')