% This file should run after the main gendata* file. It will generate the 
% SHELFICEloadAnomalyFile = 'phi0surf.exp1.jmd95z',
% in data.shelfice

% flag for construct g.HFacC, or use grid_glued.nc
% there is a small difference between the constructed result and the model
% grid result, in the bottom corners of the ice shelf
% probably again related to partial cells
use_model_grid=0; % if 0, construct custom pressure load file

if local_plot
    use_model_grid=0;
end

% equation of state
%eos = 'linear';
eos = 'jmd95z';
%eos = 'mdjwf';

% custom grid, salinity, temp

zgp1 = [0,cumsum(dzfile)];
zc = .5*(zgp1(1:end-1)+zgp1(2:end)); % center of cell
zg = zgp1(1:end-1); % top of cell
zb=zgp1(2:end); % bottom of cell
dz = diff(zgp1);

%generate_stratification_CDW26;
tref = -1.80*ones(nz,1);

% Gravity
gravity=9.81;
rhoConst = 1030;

% compute potential field underneath ice shelf
talpha = 2e-4;
sbeta  = 7.4e-4;
t    = tref;
k=1;

sref = salt_range';
s    = sref;
dzm = abs([zg(1)-zc(1) .5*diff(zc)]);
dzp = abs([.5*diff(zc) zc(end)-zg(end)]);

p = abs(zc)*gravity*rhoConst*1e-4; % pressure field
dp = p;
kp = 0;

while rms(dp) > 1e-13
    phiHydF(k) = 0;
    p0 = p;
    kp = kp+1
    for k = 1:nz
        switch eos
            case 'linear'
                drho = rhoConst*(1-talpha*(t(k)-tref(k))+sbeta*(s(k)-sref(k)))-rhoConst;
            case 'jmd95z'
                drho = densjmd95(s(k),t(k),p(k))-rhoConst;
            case 'mdjwf'
                drho = densmdjwf(s(k),t(k),p(k))-rhoConst;
            otherwise
                error(sprintf('unknown EOS: %s',eos))
        end
        phiHydC(k)   = phiHydF(k) + dzm(k)*gravity*drho/rhoConst;
        phiHydF(k+1) = phiHydC(k) + dzp(k)*gravity*drho/rhoConst;
    end
    switch eos
        case 'mdjwf'
            p = (gravity*rhoConst*abs(zc) + phiHydC*rhoConst)/gravity/rhoConst;
    end
    dp = p-p0;
end


if use_model_grid==0
    % construct approximation of g.HFacC, only used to get msk
    g.HFacC=zeros(length(dxfile),length(dyfile),length(dzfile));

    for i=1:nz
        temp_HFacC_array=zeros(size(bathy_combined));
        temp_HFacC_array(abs(zb(i))>abs(icetopo) & abs(zg(i))<abs(bathy_combined))=1;
        % bottom of grid cell is deeper than icetopo
        % top of grid cell is shallower than bathy
        % ==> water in that cell ==> assign value 1

            % possibly an error here related to partial cells
        g.HFacC(:,:,i)=temp_HFacC_array;
    end

    % the rest is the same but the way of accounting for partial cells is
    % switched
    msk=sum(g.HFacC,3); msk(msk>0)=1;
    phi0surf = zeros(size(bathy_combined));
    for ix=1:size(bathy_combined,1)
        for iy=1:size(bathy_combined,2)
            k=find(abs(zg)<abs(icetopo(ix,iy)), 1, 'last' );
            if isempty(k)
                k=0;
            end
            if k>0
                kp1=min(k+1,nz);

                % account for partial cells
                %drloc=1-g.HFacC(ix,iy,k); % for grid_glued.nc method
                drloc=(abs(icetopo(ix,iy))-abs(zg(k)))/dz(k); % for general method
                % probably should set drloc=0 if less than threshold,
                % see docs 2.11.6

                % do integration
                dphi = phiHydF(kp1)-phiHydF(k);
                phi0surf(ix,iy) = (phiHydF(k)+drloc*dphi)*rhoConst*msk(ix,iy);
            end
        end
    end
end

if use_model_grid==1
    g.HFacC=ncread('grid_glued.nc','HFacC');
    msk=sum(g.HFacC,3); msk(msk>0)=1;
    phi0surf = zeros(size(bathy_combined));
    for ix=1:size(bathy_combined,1)
        for iy=1:size(bathy_combined,2)
            k=find(abs(zg)<abs(icetopo(ix,iy)), 1, 'last' );
            if isempty(k)
                k=0;
            end
            if k>0
                kp1=min(k+1,nz);
                drloc=1-g.HFacC(ix,iy,k);
                %drloc=(abs(icetopo(ix,iy))-abs(zg(k)))/dz(k);
                dphi = phiHydF(kp1)-phiHydF(k);
                phi0surf(ix,iy) = (phiHydF(k)+drloc*dphi)*rhoConst*msk(ix,iy);
            end
        end
    end
end

figure;
h=pcolor(phi0surf);
set(h, 'EdgeColor', 'none');
title('Shelf ice pressure load file')
colorbar

% write file
fid=fopen(['phi0surf.exp1.' eos],'w','b'); fwrite(fid,phi0surf,acc);fclose(fid);

