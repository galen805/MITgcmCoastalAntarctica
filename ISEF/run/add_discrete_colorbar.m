% adds discrete colorbar to a tiled layout
function []=add_discrete_colorbar(t1,width,levels,cmap,label,tick_precision,tick_frequency)

% first add the colormap to the main axes
colormap(gca, cmap);
clim([levels(1),levels(end)])

% then go to colorbar axes
cb_ax = nexttile(t1,width);
axes(cb_ax);


numLevels=length(levels);
[Yc, Xc] = meshgrid(0:numLevels+1, [0 1]);
Z = repmat(0:numLevels+1, 2, 1);
% pcolor(cb_ax, Xc, Yc, Z,'linewidth',1,'edgecolor','k');
h = pcolor(cb_ax, Xc, Yc, Z);
set(h, 'EdgeColor', 'k', 'LineWidth', 1);
set(cb_ax, 'YDir', 'normal', ...
           'YAxisLocation', 'right', ...
           'XTick', [], ...
           'YTick', 1:tick_frequency:numLevels, ...
           'YTickLabel', string(round(levels(1:tick_frequency:end),tick_precision)), ...
           'Box', 'on');
ylim(cb_ax, [1 numLevels+1]);
pbaspect(cb_ax, [1 15 1]);
ylabel(cb_ax, label);
colormap(cb_ax, cmap);
set(cb_ax,'fontsize',20)