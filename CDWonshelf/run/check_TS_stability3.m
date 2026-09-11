% check stratification to ensure stable water column
%generate_stratification_CDW26

% using same EOS as MITgcm
eos = 'jmd95z'; % unused, just for reference
% shelf ice bin file has switch for this

gravity=9.81;
rhoConst = 1030;

%p = abs(zc)*gravity*rhoConst*1e-4; % pressure field
p=0.1*ones(size(zc)); % use surface pressure; gives potential density which is needed for N2
%s=salt_range';

% get 2D N2 field at center of domain
N2=zeros(length(latc),length(zc)-1);
slice=nx/2;
for j=1:length(latc)
    t=init_T_file(slice,j,:);
    s=init_S_file(slice,j,:);
    for k=1:length(zc)-1
        if -1*zc(k)<bathy_combined(slice,j)
            N2(j,k)=NaN;
        elseif -1*zc(k)>icetopo(slice,j)
            N2(j,k)=NaN;
        else
            rhok=densjmd95(s(k),t(k),p(k));
            rhokp1=densjmd95(s(k+1),t(k+1),p(k+1));
            dz=-(zc(k)-zc(k+1));
            drhodz=(rhok-rhokp1)/dz;
            N2(j,k)=-gravity/rhoConst*drhodz;
    
        end
    end
end

minimum_stability=min(N2,[],'all')

figure;
h=pcolor(latc,-zc(2:end),N2');
set(h, 'EdgeColor', 'none');
cb=colorbar;
title('Stability Nsq')
%saveas(gcf,'N2_stab.png')