function blasius_solver()
    % Define the Blasius equation as a system of first-order ODEs
    blasius_eq = @(eta, y) [y(2); y(3); -0.5 * y(1) * y(3)];

    % Set the boundary conditions
    eta0 = 0;
    y0 = [0; 0; 0.332];

    % Define the integration interval
    eta_span = [eta0 10];

    % Solve the ODE using ode45
    [eta, solution] = ode45(blasius_eq, eta_span, y0);

    % Extract the solution components
    theta = solution(:, 1);
    u = solution(:, 2);
    f = solution(:, 3);

    % Plot the results
    figure;
    % subplot(3, 1, 1);
    % plot(theta, eta);
    % title('Blasius Equation - Boundary Layer Thickness (\theta)');
    % ylabel('\eta');

    % subplot(3, 1, 2);
    plot(u, eta);
    title('Blasius Equation - Velocity Profile (u)');
    ylabel('\eta');
    xlabel('u/U_e')

    % subplot(3, 1, 3);
    % plot(f, eta);
    % title('Blasius Equation - Shear Stress Profile (f)');
    %ylabel('\eta');

end
