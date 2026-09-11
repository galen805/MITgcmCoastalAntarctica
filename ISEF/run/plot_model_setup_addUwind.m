%% plot model setup

% wind in x-y view
figure
m_proj('Transverse mercator','lon',[min(lonc) max(lonc)],'lat',[min(latc) max(latc)]);

% vwind1(bathy_combined==0) = NaN;
vwind_to_plot = vwind(:,:,340);
vwind_to_plot(bathy_combined==0) = NaN;
vwind_to_plot(icetopo<0) = NaN;

m_pcolor(lonc,latc,vwind_to_plot');
colormap('jet')

shading flat;
colorbar
caxis([0 wind_main]);

hold on
m_grid('box','fancy','tickdir','out');

xlabel('Longitude ($^{\circ}E$)','interpreter','latex');
ylabel('Latitude ($^{\circ}N$)','interpreter','latex');
title(['V wind field (m/s)'],'interpreter','latex','Fontsize',20);
% print('-depsc',['S_overview_z_',num2str(Z(depth_in)),'_',num2str(exact_time),'days.eps']);
saveas(gcf,['wind_field.png']);


% ice shelf topography in y-z view
figure
plot(latc, ice_line, 'linewidth',2);
xlabel('latitude','fontsize',16)
ylabel('depth (m)','fontsize',16)
title(['ice shelf shape in y-z view'],'interpreter','latex','Fontsize',20);
axis tight

saveas(gcf,['iceshelftopo_side.png']);


% iceshelf topography in x-y view
figure
m_proj('Transverse mercator','lon',[min(lonc) max(lonc)],'lat',[min(latc) max(latc)]);
% m_proj('lambert','lon',[min(lonc) max(lonc)],'lat',[min(latc) max(latc)]);

icetopo_plot2d = icetopo;
icetopo_plot2d(icetopo_plot2d==0) = NaN;

m_pcolor(lonc,latc,icetopo_plot2d');
%caxis([-400 0])
colormap('jet')

shading flat;
colorbar

hold on
m_grid('box','fancy','tickdir','out');

xlabel('Longitude ($^{\circ}E$)','interpreter','latex');
ylabel('Latitude ($^{\circ}N$)','interpreter','latex');
title(['ice shelf shape (m)'],'interpreter','latex','Fontsize',20);
saveas(gcf,['iceshelftopo.png']);


% U-wind in x-y view
figure
m_proj('Transverse mercator','lon',[min(lonc) max(lonc)],'lat',[min(latc) max(latc)]);
% 
% uwind1(bathy_combined==0) = NaN;
% uwind2(bathy_combined==0) = NaN;
% %uwind(bathy_combined==0) = NaN;
uwind_to_plot = uwind(:,:,340);
uwind_to_plot(bathy_combined==0) = NaN;
uwind_to_plot(icetopo<0) = NaN;

m_pcolor(lonc,latc,uwind_to_plot');
%m_pcolor(lonc,latc,uwind2');
%m_pcolor(lonc,latc,uwind(:,:,605)')
colormap('jet')

shading flat;
colorbar
% caxis([0 20]);

hold on
m_grid('box','fancy','tickdir','out');

xlabel('Longitude ($^{\circ}E$)','interpreter','latex');
ylabel('Latitude ($^{\circ}N$)','interpreter','latex');
title(['U wind field (m/s)'],'interpreter','latex','Fontsize',20);
% print('-depsc',['S_overview_z_',num2str(Z(depth_in)),'_',num2str(exact_time),'days.eps']);
saveas(gcf,['U_wind_field.png']);