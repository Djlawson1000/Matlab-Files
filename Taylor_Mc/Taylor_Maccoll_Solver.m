%% Taylor-Maccoll Solver
% Following the John Anderson numerical solution strategy for solving the Taylor-Maccoll 3D flow over a sharp cone
% This approach invloves an inverse method, where the shock angle is chosen
% and the cone that supports that shock is calculated

% Free Stream Mach Number, M1

M1 = 5;

% Assumed Shock Wave Angle in degrees (converted to radians below), B

B = deg2rad(30);

% Ratio of Specifc Heats, G

G = 1.4;

%% Step 1) Using Oblique shock relations from Anderson section 4.13 to calculate Mach Number just past the shock (M2) and the deflection angle (D)

D = atan((2*cot(B))*((M1.^2)*(sin(B).^2)-1)./(M1.^2*(G + cos(2*B)) + 2));

M1n = M1*sin(B);

M2n = ((2 + (G -1)*M1n.^2)./(2*G*(M1n.^2) - (G - 1))).^0.5;

M2 = M2n./sin(B - D);

%% Step 2) Using M2 and D, find the normal non-dimensional velocity components Vr (radial) and Vt (theta) from the geometry of Fig 10.4 in Anderson and V (Eq. 10.16)

V = ((2./((G - 1)*(M2.^2))) + 1).^-0.5;

Vr = V*(cos(B - D));

Vt = -V*(sin(B - D));

%% Step 3) Using an ODE23 function, calculate the required Cone angle (theta_C) to support the M1 and B conditions

list = fliplr(linspace(0.001, B, 100));

options = odeset('Events', 'on');
[theta,v] = ode23('taylor_maccoll_eqn', list, [Vr, Vt], options);

theta_C_rad = theta(end);

Theta_C_Deg = rad2deg(theta(end))

Vr_Cone = V*(cos(theta_C_rad))





