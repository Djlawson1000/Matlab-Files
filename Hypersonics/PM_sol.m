function [P_2, T_2, M2] = PM_sol(theta, M1, P, T, gamma)
    if nargin < 5
        gamma = 1.362;
    end
    gamma = 1.362;
    M2 = get_M2_prandtl_meyer(M1, deg2rad(abs(theta)));
    T_2 = T * (1 + 0.5 * (gamma - 1) * M1^2) / (1 + 0.5 * (gamma - 1) * M2^2);
    P_2 = P * ((1 + (gamma - 1) / 2 * M1^2) / (1 + (gamma - 1) / 2 * M2^2))^(gamma / (gamma - 1));
end

