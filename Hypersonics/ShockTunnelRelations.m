%% Description: Shock Tunnel Relations
% Nicholas Hawley
% Adapted from 'ShockTubeRelations.m' by Timothy Moren and 
% 'ShockTubesRelation_Matt.mat' by Matt Orcutt. 
%
% Note 1: Equation reference numbers conform to those in Modern Compressible 
% Flow (3rd Edition) by Anderson
% 
% Note 2: A nozzle with a double-ramp test model is added after shock tube
% relations. The double-ramp model includes external matlab functions
% 'obliqueShock.m' 'PMexpansion_M2out.m' 'waveAngle.m' and to perform those
% calculations. They are not needed for the nozzle expansion calculations
% 
% Note 3: 
% 'obliqueShock.m' 'PMexpansion_M2out.m' 'waveAngle.m' and to perform those
% 
% When comparing to the WiSTL code, the relevant parameters are P1, T1, T4, 
% and whatever value of Ms is calculated. 
% 2/16/2024
close all; clear; clc;

%% Initialization: Driver, Driven, and Nozzle Conditions

% Driver Tube (station 4), Gas = Air
p4 = 82; % [psi]
t4 = 293; % [K]
% t4 = 4.775944444444444e+02;
gamma4 = 1.4;
r4 = 287; % [J/(kg*K)]
a4 = sqrt(gamma4 * r4 * t4); % [m/s]

% % Driver Tube (station 4), Gas = He
% p4 = 250; % [psi]
% t4 = 293; % [K]
% gamma4 = 5/3;
% r4 = 2077.1; % [J/(kg*K)]
% a4 = sqrt(gamma4 * r4 * t4); % [m/s]

% Driven Tube (station 1), Gas = Air
p1 = 14.7; % [psi]
% air_temp_1 = 400; % [F], this can be heated or room temp
% t1 = (air_temp_1 - 32) * (5/9) + 273.15; % [K]
t1 = 293;
gamma1 = 1.4;
r1 = 287; % [J/(kg*K)]
a1 = sqrt(gamma1 * r1 * t1); % [m/s]
l1 = 5; % [m]

% CPG assumption -- note that this assumption applies for stations 3/4 and
% 1/2. I.e., use gamma1 for calculating ratios across stations 1 and 2,
% and use gamma4 for calculating ratios across stations 3 and 4.

% Ratios
p4p1 = p4 / p1;
t4t1 = t4 / t1;
a4a1 = a4 / a1;
gg_1 = (gamma1 + 1) / (gamma1 - 1);
g2g_1 = (2 * gamma1) / (gamma1 + 1);
gg_4 = (gamma4 + 1) / (gamma4 - 1);
g2g_4 = (2 * gamma4) / (gamma4 + 1);

% Nozzle Design Mach Number
Mdesign = 4;


%% Calculations: Incident Shock
% find pressure and temperature ratios across the incident shock
    % solve for p2/p1 using eq. 7.94
    p2p1_func = @(p2p1_ratio) p2p1_ratio*((1 - ((gamma4 - 1) * (1/a4a1) * ...
        (p2p1_ratio - 1)) / (sqrt(2 * gamma1 * (2 * gamma1 + (gamma1 + ...
        1) * (p2p1_ratio - 1)))))^(-2 * gamma4 / (gamma4 - 1))) - p4p1;
    p2p1_guess = 1;
    p2p1 = fsolve(p2p1_func, p2p1_guess, optimset('Display', 'off'));

    % solve for t2/t1 using eq. 7.10
    t2t1 = p2p1 * ((gg_1 + p2p1) / (1 + gg_1 * p2p1));

% find wave velocity, W, of the incident shock
Wi = a1 * sqrt((1 / g2g_1) * (p2p1 - 1) + 1); % [m/s]

% find Mach number of the incident shock, Ms using eq. b/t 7.11 and 7.12
Ms = Wi / a1; % gives the same result as eq. 7.13 below
% Ms = sqrt(g2g * (p2p1 - 1) + 1); % eq. 7.13

% find induced velocity up behind the wave using eq. 7.16
up = (a1 / gamma1) * (p2p1 - 1) * ((g2g_1 / (p2p1 + (1/gg_1)))^0.5); %[m/s]

% find induced particle mach number using eq. 7.17
Mp = (up / a1) * sqrt(1 / t2t1);

% Check p4/p1 with eq. 8 from GALCIT-6 CalTech
check = (1 + g2g_1 * (Ms^2 - 1)) / ((1 - (gamma4 - 1) / (gamma1 + 1) * ...
    a1 / a4 * (Ms^2 - 1) / Ms)^(2 * gamma4 / (gamma4 - 1)));

%% Calculations: Reflected Shock 
% find reflected shock Mach number using eq. 7.23
Mr_func = @ (Mr_var) (Ms / (Ms^2 - 1)) * sqrt(1 + (2 * (gamma1 - 1) / ...
    (gamma1 + 1)^2) * (Ms^2 - 1) * (gamma1 + (1 / Ms^2))) - ...
    (Mr_var / (Mr_var^2 - 1));
Mr = fzero(Mr_func, Ms, optimset('Display', 'off'));

% calculate pressure ratio after reflected shock, p5/p2 adapting eq. 7.12
% using the reflected shock Mach (Mr) instead of incident shock Mach (Ms)
    % Note: refere to eq. 3.57 where M1 = Mr
p5p2 = 1 + g2g_1 * (Mr^2 - 1);
% check = g2g_1 * ((2 * gamma1 * Ms^2 + 1 - gamma1) / ...
%     ((gamma1 - 1) * Ms^2 + 2)) - (1 / gg_1); % eq. 16 from GALCIT-6 CalTech

% calculate pressure ratio p5/p1 from eq. 10 of GALCIT-6 CalTech
p5p1 = 1 + 2 * (p2p1 - 1) * ((1 + (0.5 + 1 / gg_1) * (Ms^2 - 1)) / ...
    (1 + 1 / gg_1 * (Ms^2 - 1)));
% check = p5p2 * p2p1;

% calculate temperature ratio, T5/T2 adapting eq. 3.59 where M1 = Mr
t5t2 = (1 + g2g_1 * (Mr^2 - 1)) * ((2 + (gamma1 - 1) * Mr^2) / ...
    ((gamma1 + 1) * Mr^2));

% Calculate speed of sound ahead of the reflected shock
a2 = sqrt(gamma1 * r1 * t2t1 * t1); % [m/s]

% Calculate wave speed of the reflected shock: "velocity of the gas behind
% the shock wave relative to the wave
    % note that Wr + up_shock is the velocity of the gas ahead of the 
    % reflected shock wave relative to the wave. 
    % SEE PAGE 274 OF ANDERSON FOR A PICTURE AND BOTTOM OF 275 FOR EQUATION
Wr = Mr * a2 - up; % [m/s]

% calculated reflected shock wave speed using eq. 11 from GALCIT-6 CalTech
    % note: can calculate using eq. ____ in anderson
ur = a1 * ((1 + 2 * (1 / gg_1) * (Ms^2 - 1)) / Ms); % [m/s]

%% Calculations: Incident Expansion Wave 
% Section 7.7 in Anderson third edition
% calculate p3/p4 using pressure ratios (p3 = p2), p. 298
    % Note: mulitple ways to calculate p3/p4 -- could use eq. 7.90 where
    % u3 = up_shock
p3p4 = p2p1 / p4p1;
% p3p4 = (1 - ((gamma4 - 1) / 2) * (up_shock / a4))^(2 * gamma4 / ...
    % (gamma4 - 1));

% calculate rho3/rho4 using isentropic relations (p. 298 anderson)
rho3rho4 = p3p4^(1 / gamma4);

% calcualte t3/t4 using isentropic relations (p. 298 anderson)
t3t4 = p3p4^((gamma4 - 1)  / gamma4);


%% Calculations: Reflected Expansion Wave

% initialize 
    % Note: all equations found in sections 7.6 and 7.7 of anderson
syms x_new
nwaves = 5; % number of expansion waves
a = zeros(nwaves); % sound speed of wave
u = zeros(nwaves); % velocity of wave; notes, velocity at wall is 0
u(1,1:nwaves) = linspace(0,up,nwaves); % first row is the u velocity of the C+ line
jplus = zeros(nwaves); % right running riemann invariant
jminus = zeros(nwaves); % left running riemann invariant
x = zeros(nwaves,nwaves+1); % x-coord
t = zeros(nwaves,nwaves+1); % t vals

% first point -- point at which first wave hits wall and is reflected
a(1,1) = a4;
% u(1,1) is already zero
jplus(1,1) = u(1,1) + 2 * a(1,1) / (gamma4 - 1); % eq. 7.73
jminus(1,1) = -jplus(1,1); % eq. 7.74 (kind of)
x(1,1) = 0;
t(1,1) = -l1 / (u(1,1) - a(1,1));

% solve the points that are along the first C+ characteristic line
for j = 2:nwaves
    a(1,j) = a(1,1) - (gamma4 - 1) / 2 * u(1,j); % eq. 7.84
    jplus(1,j) = u(1,j) + 2 * a(1,j) / (gamma4 - 1); % eq. 7.73
    jminus(1,j) = u(1,j) - 2 * a(1,j) / (gamma4 - 1); % eq. 7.74

    cminus = 1 / (u(1,j) - a(1,j)); % fig. 7.13
    cplus = tan(0.5 * (atan(1 / (u(1,j-1) + a(1,j-1))) + ...
        atan(1 / (u(1,j) + a(1,j))))); % bottom eq. on p. 296 
    
    % solve for x setting up system of equations using two point-slope form
    % equations with points 1 and points 0
    x(1,j) = solve(t(1,j-1) + cplus * (x_new - x(1,j-1)) == cminus * ...
        (x_new - l1),x_new);
    t(1,j) = t(1,j-1) + cplus * (x(1,j) - x(1,j-1));

end

% solve all other points 
for i = 2:nwaves
    % points on the wall -- i.e. diagonal elements of matrix
    u(i,i) = 0; % wall boundary condition u4 = 0
    jminus(i,i) = jminus(i-1,i); % invariant is same along same C- line
    jplus(i,i) = -jminus(i,i);
    a(i,i) = (gamma4 - 1) / 4 * (jplus(i,i) - jminus(i,i)); % eq. 7.75

    cminus = tan(0.5 * (atan(1 / (u(i-1,i) - a(i-1,i))) + ...
        atan(1 / (u(i,i) - a(i,i))))); % top eq. on p. 296
    x(i,i) = 0;
    t(i,i) = t(i-1,i) + cminus * (x(i,i) - x(i-1,i));
    
    % every other point
    for j = i+1:nwaves
        jminus(i,j) = jminus(i-1,j);
        jplus(i,j) = jplus(i,j-1);
        u(i,j) = (jplus(i,j) + jminus(i,j)) / 2;
        a(i,j) = (gamma4 - 1) / 4 * (jplus(i,j) - jminus(i,j));

        cplus = tan(0.5 * (atan(1 / (u(i,j-1) + a(i,j-1))) + ...
            atan(1 / (u(i,j) + a(i,j)))));
        cminus = tan(0.5 * (atan(1 / (u(i-1,j) - a(i-1,j))) + ...
            atan(1 / (u(i,j) - a(i,j)))));

        x(i,j) = solve(t(i-1,j) + cminus * (x_new - x(i-1,j)) == ...
            t(i,j-1) + cplus * (x_new - x(i,j-1)),x_new); % point slope sys of eq
        t(i,j) = t(i,j-1) + cplus * (x(i,j) - x(i,j-1)); % point slope sys of eq

    end
end

% shift wall to -l1
x = x - l1;

% Plot points of intersection of each reflected wave
figure(3)
plot(x,t,'ko');
xlabel('x [m]');
ylabel('t [sec]');

% plot wave lines

close all;

%% Calculations: x-t diagram

% TO BE COMPLETED LATER


%% Calculations: Nozzle Expansion


% static pressure after reflected shock
p5 = p5p1 * p1; % [psi]

% static temperature and speed of sound after reflected shock
t5 = t5t2 * t2t1 * t1; % [K]
a5 = sqrt(gamma1 * r1 * t5); % [m/s]

% Output station 5 conditions
disp('STATION 5 CONDITIONS (after reflected shock):');
disp(['Pressure (P5): ',num2str(p5),' [psi]']);
disp(['Temperature (T5): ',num2str(t5),' [K]',newline]);

% Now, using t5 and p5 as the stagnation (total) conditions, find static
% pressure and temperature p6 and t6 after nozzle expansion
p6 = p5 * (1 + (gamma1 - 1) / 2 * Mdesign^2)^(-gamma1 / (gamma1 -1)); %[psi]
t6 = t5 * (1 + (gamma1 - 1) / 2 * Mdesign^2)^-1; % [K]

% % % Test Model Conditions -- double ramp with reflected shock
% % Note: comment this section out if not using double ramp with reflected
% % cowl. Nick Hawley added for calculating conditions specific to his test
% % section
% % clear p6 t6 Mdesign
% % 
% % generate potential nozzle mach numbers
% % Mdesign = 1.5:0.5:8;
% % Mdesign = 2.5;
% % 
% % allocate memory
% % alpha_best = zeros(length(Mdesign),4);
% % mach_best = alpha_best;
% % pressure_best = mach_best;
% % temp_best = pressure_best;
% % beta_best = temp_best;
% % 
% % calculate conditions at each stage of test section for each mach
% % number using external function 'doubleRampOptimization.m'
% % for i = 1:length(Mdesign)    
% %     [alpha_best(i,:),mach_best(i,:),pressure_best(i,:),temp_best(i,:), ...
% %         beta_best(i,:)] = doubleRampOptimization(Mdesign(i),p5,t5,gamma1);
% % end
% % 
% % find maximum temperature in the channel (last column)
% % [temp_max,temp_max_index] = max(temp_best(:));
% % [temp_max_row, temp_max_col] = ind2sub(size(temp_best), temp_max_index);
% % disp('MAXIMUM CHANNEL TEMPERATURE AND NOZZLE MACH NUMBER:')
% % disp(['Maximum static temperature achieved in channel: ', ...
% %     num2str(temp_max),' [K]'])
% % 
% % find the optimal nozzle mach number
% % mach_nozzle = Mdesign(temp_max_row);
% % disp(['Mach number of nozzle producing highest Tchannel: ', ...
% %     num2str(mach_nozzle),' [-]']);
% % 
% % pull the corresponding geometry and test model conditions
% % alpha_model = alpha_best(temp_max_row,:); % [deg]
% % mach_model = mach_best(temp_max_row,:);
% % pressure_model = pressure_best(temp_max_row,:); % [psi]
% % temp_model = temp_best(temp_max_row,:); % [K]
% % beta_model = beta_best(temp_max_row,:); % [deg]
% % 
% % BUILD TABLE FOR CLEAN OUTPUT INTO COMMAND WINDOW
% % disp([newline,'OPTIMAL DOUBLE-RAMP GEOMETRY AND CONDITIONS:']);
% % Define the parameter names
% % parameters = {'Ramp Angle (alpha)', 'Mach Number', 'Pressure', 'Temperature', ...
% %     'Shock Angle (beta)'};
% % Define the values for each case
% % values = [alpha_model; mach_model; pressure_model; temp_model; beta_model];
% % units = {'deg', '', 'psi', 'K', 'deg'};
% % cases = {'Free Stream', 'Ramp 1', 'Ramp 2', 'Channel'};
% % max_parameter_length = max(cellfun(@length, parameters));
% % max_unit_length = max(cellfun(@length, units));
% % mach_model_cell = num2cell(mach_model);
% % max_mach_length = max(cellfun(@(x) numel(sprintf('%.2f', x)), mach_model_cell));
% % Display the header
% % fprintf('%-*s', max_parameter_length, 'Parameter');
% % for i = 1:numel(cases)
% %     fprintf('\t%-*s', max_mach_length + max_unit_length + 6, cases{i}); % Adding 6 for proper alignment with units
% % end
% % fprintf('\n');
% % Display the lines
% % line_length = max_parameter_length + 8 + (numel(cases) * (max_mach_length + max_unit_length + 6)); % 8 accounts for tabs
% % fprintf(repmat('-', 1, line_length));
% % fprintf('\n');
% % Display the data
% % for i = 1:numel(parameters)
% %     fprintf('%-*s', max_parameter_length, parameters{i});
% %     for j = 1:numel(cases)
% %         fprintf('\t%*.2f %s', max_mach_length, values(i, j), units{i});
% %     end
% %     fprintf('\n');
% % end
% % 
% % % Save data for Air, 250psi (COMMENT OUT IF NOT USING AIR AT 250 PSI)
% % % Define the data to be saved
% %     % nozzle mach (col 1), channel pressure/temp (col 4)
% % data = [mach_best(:,1) pressure_best(:,4) temp_best(:,4)];
% % 
% % % Define the headers for each column
% % headers = {'Mach', 'Pressure', 'Temperature',};
% % 
% % % convert to table
% % data_table = array2table(data,'VariableNames',headers);
% % 
% % % Write data and headers to a CSV file
% % filename = 'nozzle_conditions_air_250.csv';
% % writetable(data_table, filename); % Write data to CSV file
% % 
% % % Save data for Air, 500psi (COMMENT OUT IF NOT USING AIR AT 500 PSI)
% % % Define the data to be saved
% %     % nozzle mach (col 1), channel pressure/temp (col 4)
% % data = [mach_best(:,1) pressure_best(:,4) temp_best(:,4)];
% % 
% % % Define the headers for each column
% % headers = {'Mach', 'Pressure', 'Temperature',};
% % 
% % % convert to table
% % data_table = array2table(data,'VariableNames',headers);
% % 
% % % Write data and headers to a CSV file
% % filename = 'nozzle_conditions_air_500.csv';
% % writetable(data_table, filename); % Write data to CSV file
% % 
% % % Save data for He, 250psi (COMMENT OUT IF NOT USING He AT 250 PSI)
% % % Define the data to be saved
% %     % nozzle mach (col 1), channel pressure/temp (col 4)
% % data = [mach_best(:,1) pressure_best(:,4) temp_best(:,4)];
% % 
% % % Define the headers for each column
% % headers = {'Mach', 'Pressure', 'Temperature'};
% % 
% % % convert to table
% % data_table = array2table(data,'VariableNames',headers);
% % 
% % % Write data and headers to a CSV file
% % filename = 'nozzle_conditions_he_250.csv';
% % writetable(data_table, filename); % Write data to CSV file
% % 
% % % Save data for He, 500 psi (COMMENT OUT IF NOT USING He AT 500 PSI)
% % % Define the data to be saved
% %     % nozzle mach (col 1), channel pressure/temp (col 4)
% % data = [mach_best(:,1) pressure_best(:,4) temp_best(:,4)];
% % 
% % % Define the headers for each column
% % headers = {'Mach', 'Pressure', 'Temperature'};
% % 
% % % convert to table
% % data_table = array2table(data,'VariableNames',headers);
% % 
% % % Write data and headers to a CSV file
% % filename = 'nozzle_conditions_he_500.csv';
% % writetable(data_table, filename); % Write data to CSV file
% % 
% % % Save data for Air, 250 psi, Vacuumed (COMMENT OUT IF NOT USING)
% % % Define the data to be saved
% %     % nozzle mach (col 1), channel pressure/temp (col 4)
% % data = [mach_best(:,1) pressure_best(:,4) temp_best(:,4)];
% % 
% % % Define the headers for each column
% % headers = {'Mach', 'Pressure', 'Temperature'};
% % 
% % % convert to table
% % data_table = array2table(data,'VariableNames',headers);
% % 
% % % Write data and headers to a CSV file
% % filename = 'nozzle_conditions_air_250_vac.csv';
% % writetable(data_table, filename); % Write data to CSV file
% % 
% % % Save data for He, 250 psi, Vacuumed (COMMENT OUT IF NOT USING)
% % % Define the data to be saved
% %     % nozzle mach (col 1), channel pressure/temp (col 4)
% % data = [mach_best(:,1) pressure_best(:,4) temp_best(:,4)];
% % 
% % % Define the headers for each column
% % headers = {'Mach', 'Pressure', 'Temperature'};
% % 
% % % convert to table
% % data_table = array2table(data,'VariableNames',headers);
% % 
% % % Write data and headers to a CSV file
% % filename = 'nozzle_conditions_he_250_vac.csv';
% % writetable(data_table, filename); % Write data to CSV file
% % 
% % 
% % 
% % 
