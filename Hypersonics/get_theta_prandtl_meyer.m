function theta = get_theta_prandtl_meyer(M1, M2, gamma)
    gamma = 1.4;
    vm1 = sqrt((gamma + 1) / (gamma - 1)) * atan(sqrt((gamma - 1) / (gamma + 1) * (M1^2 - 1))) - atan(sqrt(M1^2 - 1));
    vm2 = sqrt((gamma + 1) / (gamma - 1)) * atan(sqrt((gamma - 1) / (gamma + 1) * (M2^2 - 1))) - atan(sqrt(M2^2 - 1));
    theta = vm2 - vm1;
    theta = rad2deg(theta);
end
