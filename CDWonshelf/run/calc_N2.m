function N2_range = calc_N2(dens_range,z)
% calculate 1D stratification


gravity=9.81;
rhoConst = 1030;

drhodz = gradient(dens_range, -1*abs(z));  % dρ/dz
N2_range = -gravity / rhoConst * drhodz;