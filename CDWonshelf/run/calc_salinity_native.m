function S_range = calc_salinity_native(dens_range,T_range)
% calculate salinity from 1D temp and density using native MITgcm eos

nz=length(dens_range);
% calculate salinity profile
p=0.1; % surface pressure
for k=1:nz
    n_iter=0;
    rescheck=1;

    match_dens=dens_range(k)+1000;

    % Initial guess
    S = 33.0;

    while rescheck > 1e-5
        % Compute current density and residual
        dens = densjmd95(S, T_range(k), p);
        res = dens - match_dens;

        % Finite difference approximation of d(dens)/dS
        dS = 1e-4;
        dens_dS = densjmd95(S + dS, T_range(k), p);
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
            fprintf('Too many iterations at level %d\n', k);
            break;
        end
    end

    S_range(k) = S;

end

end