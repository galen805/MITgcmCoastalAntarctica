function add_coast_func(X,Y,bathy_combined)
% plots coastline geometry
% load plot_model_data.mat % call in main script

% Time: 2023.06.11
% Use m_patch to indicate the coast.
% Previous version uses m_pcolor,
% but it might not work well with system macOS 12 and later versions.
% When freezeColors function is used after m_pcolor,
% the code will freeze in macOS 12 and macOS 13.

mask_depth = bathy_combined; % use bathy to get the boundaries of the water area
    % need to update this for when bathy intersects ice shelf

mask_depth_boundary_cell = bwboundaries(mask_depth==0);
mask_depth_boundary = mask_depth_boundary_cell{1};
mask_depth_boundary_column_1 = mask_depth_boundary(:,1);
mask_depth_boundary_column_2 = mask_depth_boundary(:,2);

% m_pcolor(X, Y, mask_depth'); shading flat;
%m_patch(X(mask_depth_boundary_column_1),Y(mask_depth_boundary_column_2),[170 170 170]/255, 'EdgeColor',[170 170 170]/255);
m_patch(X(mask_depth_boundary_column_1),Y(mask_depth_boundary_column_2),[220,220,220]/255, 'EdgeColor',[0,0,0]/255,'linewidth',2);
hold on;

% flags
plot_ice_tongue = 0;
plot_with_cape = 0;

if plot_ice_tongue == 1
    plot_ice_tongue_100=1;
    plot_ice_tongue_50=~plot_ice_tongue_100;

    if plot_ice_tongue_100 == 1
        Y_north_in  = 200;
    end
    if plot_ice_tongue_50 == 1
        Y_north_in = 150;
    end
    m_line([X(227) X(227)],[Y(101) Y(Y_north_in)],'linewi',2,'color','k','linestyle','--'); hold on;
    m_line([X(227) X(252)],[Y(Y_north_in) Y(Y_north_in)],'linewi',2,'color','k','linestyle','--'); hold on;
    m_line([X(252) X(252)],[Y(Y_north_in) Y(101)],'linewi',2,'color','k','linestyle','--'); hold on;
    m_line([X(227) X(252)],[Y(101) Y(101)],'linewi',2,'color','k','linestyle','--'); hold on;
end

% 20210407 - add headland similar to Cape Washington
if plot_with_cape % B5 case with B2 turned on
    add_cape_lines;

    % label
    % m_text(4.5,-75.6,'Headland','color','k','fontsize',28);
    %     m_text(4.5,-75.6,'Headland','color','k','fontsize',14)
end

