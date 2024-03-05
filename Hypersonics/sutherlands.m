function mu = sutherlands(T, mu_ref, Tref, S)
    if nargin < 4
        S = 110.4; % Default value for S
    end
    if nargin < 3
        Tref = 273.15; % Default value for Tref
    end
    if nargin < 2
        mu_ref = 1.716e-5; % Default value for mu_ref
    end
    
    mu = mu_ref * ((Tref + S) / (T + S)) * ((T / Tref) ^ 1.5);
end
