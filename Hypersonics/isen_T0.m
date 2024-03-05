function T0 = isen_T0(T_inf, M, gamma)
    if nargin < 3
        gamma = 1.4; % Default value for gamma, adjust if needed
    end
    
    T0 = (1 + M^2 * (gamma - 1) / 2) * T_inf;
end
