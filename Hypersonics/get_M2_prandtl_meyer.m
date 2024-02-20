function M2 = get_M2_prandtl_meyer(M1, theta, gamma)
    % Define the function for fzero
    function result = func(M2, M1, theta, gamma)
        result = ((sqrt((gamma + 1) / (gamma - 1)) * atan(sqrt((gamma - 1) / (gamma + 1) * (M2^2 - 1))) - atan(sqrt(M2^2 - 1))) - (sqrt((gamma + 1) / (gamma - 1)) * atan(sqrt((gamma - 1) / (gamma + 1) * (M1^2 - 1))) - atan(sqrt(M1^2 - 1))) - theta);
    end
    gamma = 1.362;
    % Use fzero to find the root
    options = optimset('Display', 'off');
    sol = fzero(@(M2) func(M2, M1, theta, gamma), [M1, 1e3], options);

    % Return the result
    M2 = sol;
end

