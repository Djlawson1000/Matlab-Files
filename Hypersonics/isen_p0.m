function P0 = isen_p0(P_inf, M, gamma)
    if nargin < 3
        gamma = 1.4; % Default value for gamma, adjust if needed
    end
    
    P0 = (1 + M^2 * (gamma - 1) / 2)^(gamma / (gamma - 1)) * P_inf;
end
