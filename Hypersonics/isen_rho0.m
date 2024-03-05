function rho0 = isen_rho0(rho_inf, M, gamma)
    if nargin < 3
        gamma = 1.4; % Default value for gamma
    end
    
    rho0 = (1 + M^2 * (gamma - 1) / 2)^(1 / (gamma - 1)) * rho_inf;
end
