function beta = wave_angle(M, theta, gamma)
    if nargin < 3
        gamma = 1.4;
    end
    gamma = 1.4;
    % Define the function to solve
    func = @(beta) tan(deg2rad(theta)) - 2/tan(beta) * ((M * sin(beta))^2 - 1) / (M^2 * (gamma + cos(2 * beta)) + 2);

    % Start with an initial guess for beta (e.g., 5 degrees in radians)
    initial_guess = asin(1/M);

    % Solve for beta using fsolve
    beta = fsolve(func, initial_guess);

    % Display the result
    fprintf('Computed beta: %.4f degrees\n', rad2deg(beta));
end
