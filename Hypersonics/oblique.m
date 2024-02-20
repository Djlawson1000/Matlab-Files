function [P_2, T_2, M2, beta_deg] = oblique(theta, M1, P, T, gamma)
    if nargin < 5
        gamma = 1.362;
    end
    gamma = 1.362;
    beta = wave_angle(M1, theta);
    % beta_rad = deg2rad(beta);
    M1_n = M1 * sin(beta);
    M2_n = sqrt((1 + 0.5 * (gamma - 1) * M1_n^2) / (gamma * M1_n^2 - 0.5 * (gamma - 1)));
    M2 = M2_n / sin(beta - deg2rad(theta));

    P_2 = P * (1 + gamma * M1_n^2) / (1 + gamma * M2_n^2);
    T_2 = T * (1 + 0.5 * (gamma - 1) * M1_n^2) / (1 + 0.5 * (gamma - 1) * M2_n^2);

    beta_deg = rad2deg(beta);
end
