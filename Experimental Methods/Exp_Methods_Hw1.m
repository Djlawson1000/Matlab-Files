%% Experimental Methods HW #1 

%% Problem 2 Confidence Interval

% Plot both PDF and CDF on the same graph given mean and sigma
P_bar = 20;

std = 0.5;

x = 18:0.1:22;

z = normpdf(x, P_bar, std);

plot(x,z,'-k'); hold on;

N = normcdf(x, P_bar, std); plot(x,N,'--r');

%% Probabilities of varying results

% Probability (%) measured value is less than 19.2 psi?

Prob_b = 100*normcdf(19.2, 20, 0.5);

% Probability (%) measured value is between 19.6 and 20.9 psi?

Prob_c = 100*(normcdf(20.9, 20, 0.5) - normcdf(19.6, 20, 0.5));

% Probability (%) measured value is larger than 20.5 psi?

Prob_d = 100*(1 - normcdf(20.5, 20, 0.5));

%% Problem 3 Sample Data

% Plot Meaningful histograms for the given prssure and temp data

P = [4.95, 4.87, 5.1, 4.99, 5.06, 5.09, 4.95, 5.01, 4.98, 5.01, 4.96, 4.99, 5.00, 5.05, 5.00, 5.02, 5.01, 5.06, 4.99, 5.01, 4.98];
T = [21.2, 21.3, 21.0, 21.2, 21.5, 20.9, 21.2, 21.4, 21.1, 21.2, 21.4, 21.1, 21.4, 20.9, 21.5, 21.2, 21.6, 21.1, 20.8, 20.7, 20.9];

P_h = histogram(P, 7);
T_h = histogram(T, 5);

