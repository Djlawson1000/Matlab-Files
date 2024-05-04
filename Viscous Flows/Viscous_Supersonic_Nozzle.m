

%% Pressure Profile

clear; clc; close all

%PRESSURE PROFILE__________________________________________________________
% Begin computational iteration over grid
%Pressure Distibution
x = (.0236:.0001:0.3);
P_Po = 1 - .979 * exp(-(9.7*10^(-3))./(x-0.018));
xadd = [ 0 .005 .00725 .01 .0125 .015 .0175 .02 .0225 .023 ];
xadd(end) = xadd(end)-.00034;clo
xplot= [xadd x]/100;
extraP = [ 1 .99511 .98961 .97986 .96775 .95499 .93171 .90519 .86861 .85931];
P_Po2 = [extraP, P_Po];
%Plot profile
figure
plot(xplot,P_Po2,'-','linewidth',2)

grid on
hXLabel=ylabel('$\frac{P}{P_0}$','interpreter','latex');
hYLabel=xlabel('$x$ (cm) ','interpreter','latex');
set(gca, 'FontName', 'Times')
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'LineWidth', 1, 'fontsize',20)

% Make points a better fit

% figure
% plot(xadd,extraP)
% hold on
z = polyfit(xadd,extraP,4);
xarray = linspace(min(xadd),0.024,300);
yeq = z(1).*xarray.^4 + z(2).*xarray.^3 + z(3).*xarray.^2 + z(4).*xarray.^1 + z(5);
% plot(xarray,yeq)

% Now plot the best
P_Po2 = [yeq, P_Po];
%x(1) = x(2)-.00002;
xplot= [xarray.*.97 x]/100;
%Plot profile
figure
plot(xplot,P_Po2,'-','linewidth',1)
xplot = smooth(smooth(smooth(smooth(xplot))))'; P_Po2 = smooth(smooth(smooth(smooth(P_Po2))))';
plot(smooth(smooth(xplot)),smooth(smooth(P_Po2)),'-k','linewidth',2)
grid on
hXLabel=ylabel('$\frac{P}{P_0}$','interpreter','latex');
hYLabel=xlabel('$x$ (cm) ','interpreter','latex');
set(gca, 'FontName', 'Times')
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'LineWidth', 1, 'fontsize',20)

%%
close all
%Assumption  y and psi are normal to the wall but not really true._________
% Inlet Conditions_________________________________________________________
 T_in = 2000;%K
 gamma = 1.4;
 R = 287;
 P_in = 101325;%Pa
 V_in = 340;%m/s
 a_in = sqrt(gamma*R*T_in);
 Pr_air = 0.73; % Prandtle number for air
 x_in = 0.004; % Distance across plate
 rho_in = P_in/(R*T_in);
 T_wall = 450; % K
 M_in = V_in/a_in;
rho_wall = P_in/(R*T_wall);
 
% Dynamic viscosity
 mu_wall = sqrt(T_wall).*((1.1E-21).*T_wall.^5 - (4.7E-18).*T_wall.^4+(7.6E-15).*T_wall.^3-(6.3E-12).*T_wall.^2+(3E-9).*T_wall+5.6E-7);
 nu_in = mu_wall/rho_wall;
% Cp for air function
 Cp_in = (7.9e-17*T_in^5 - 1.2e-013*T_in^4 - 3.6e-10*T_in^3 + 8.4e-7*T_in^2 - 0.00035*T_in + 1)*1000;
% K for air function
 K = 1.7e-11*T_in^3 - 5.1e-8*T_in^2 + 0.0001*T_in - 0.00053;

 % %Computational grid
j = 1;
i = 1;
x = xplot;
y_centerline = 0.0006; %m

%Create the psi space for this entire code
psi_centerline = y_centerline*V_in*rho_in;
dpsi = psi_centerline/50;
psi = (0:dpsi:psi_centerline);

% INITIAL BOUNDARY LAYER___________________________________________________
% Compute Boundary layer thickness
delta = 5.83*sqrt(nu_in*x_in/V_in);
%Relationship to blasius boundary layer
Capdelta = Pr_air^(1/3);
delta_T = delta*Capdelta;

%Initialize Variables
y_from_psi(1) = 0;         u_Uinfty(1) = 0;        eta(1) = 0;   u(1) = 0;
rho(1) = rho_wall;        TempNON_BL(1) = 1;      eta_T(1) = 0;
T(1) = T_wall;            
for i = 2:length(psi)

    if i ==2
        y_from_psi(i) = y_from_psi(i-1) + dpsi/(rho(i-1)*100);
    else
        y_from_psi(i) = y_from_psi(i-1) + dpsi/(rho(i-1)*u(i-1));
    end
    % BLASIUS BOUNDARY LAYER_______________________________________________
    % Nondim parameter to Blasius solution
    eta(i) = y_from_psi(i)./delta;
    % Blaius polynomial solution
    u_Uinfty(i) = 2*eta(i) - 2* eta(i).^3 + eta(i).^4;
    u(i) = u_Uinfty(i)*V_in;
    %TEMPERATURE BOUNDARY PROFILE__________________________________________
    % Nondim parameter to temperaure solution
    eta_T(i) = y_from_psi(i)./delta_T;
    % Nondimesnional temperature expression
    TempNON_BL(i) = 1 - 2*eta_T(i) + 2 * eta_T(i).^3 - eta_T(i).^4;
    T(i) = TempNON_BL(i)*(T_wall - T_in)+T_in;
    %Desnity change
    rho(i) = P_in/(R*T(i));

    
end

%Plot Temperature BL profile
figure
subplot(1,2,1); plot(TempNON_BL,eta_T,'linewidth',1),grid on
ylabel('\Theta','fontsize',20,'fontname','times'); xlabel('\eta_T','fontsize',20,'fontname','times')
subplot(1,2,2); plot(T,eta_T.*delta_T,'linewidth',2); grid on
ylabel('T (K)','fontsize',20,'fontname','times'); xlabel('y (m)','fontsize',20,'fontname','times')

%Plot BL profile
figure
subplot(1,2,1); plot(u_Uinfty,eta,'linewidth',1); grid on
ylabel('u/U_{\infty}','fontsize',20,'fontname','times'); xlabel('\eta','fontsize',20,'fontname','times')
subplot(1,2,2); plot(u_Uinfty.*V_in,eta.*delta,'linewidth',2); grid on
ylabel('u (m/s)','fontsize',20,'fontname','times'); xlabel('y (m)','fontsize',20,'fontname','times')

for i = 1:length(u)
    if u(i) > V_in
        u(i) = V_in;
    end
    if T(i) > T_in
        T(i) = T_in;
    end
end

%Plot Temperature BL profile
figure
plot(T,eta_T.*delta_T,'linewidth',2); grid on
xlabel('T (K)','fontsize',20,'fontname','times'); ylabel('y (m)','fontsize',20,'fontname','times')

%Plot BL profile
figure
plot(u,eta.*delta,'linewidth',2); grid on
xlabel('u (m/s)','fontsize',20,'fontname','times'); ylabel('y (m)','fontsize',20,'fontname','times')

%% Now Step through nozzle
close all
j = 1;
i = 1;

gridlengthx = length(x);
gridlengthpsi = length(psi);

%Dynamic viscosity
mu(j,:) = (1.1e-21.*T(j,:).^5 - 4.7e-18.*T(j,:).^4 + 7.6e-15.*T(j,:).^3 - 6.3e-12.*T(j,:).^2 + 3e-9.*T(j,:) + 5.6e-7).*sqrt(T(j,:)); 
%Initialize
rho(j,:) = P_in./(R.*T(j,:));
%Cp for air function
cp(j,:) = (7.9e-17*T(j,:).^5 - 1.2e-013*T(j,:).^4 - 3.6e-10*T(j,:).^3 + 8.4e-7*T(j,:).^2 - 0.00035.*T(j,:) + 1)*1000;
% K for air function
k(j,:) = 1.7e-11*T(j,:).^3 - 5.1e-8*T(j,:).^2 + 0.0001*T(j,:) - 0.00053;
%Mach number
a(j,:) = sqrt(gamma * R .* T(j,:));
M(j,:) = u(j,:)./a(j,:);
P(j,:) = ones(1,length(psi)).*P_in;

%Plot velocity and temperature profiles before iteration
figure(20)
plot(T(j,:),'linewidth',2)
hold on

% figure(21)
% plot(u(j,:),'linewidth',2)
% hold on
% ydatathis(j,:) = psi.*u(j,:).*rho(j,:);

dx = gradient(x);

P0 = 1/(P_Po2(1)/P_in);
P = P0*P_Po2;
P = [P, P(end)];
dp_dx = gradient(P(1:end-1))./dx;


for j = 1:gridlengthx
    for i = 1:gridlengthpsi
              
        if i == 1 %Wall_____________________________________________________
         
            
        elseif i <gridlengthpsi(end)%Center__________________________________
            
            A1c = dpsi; B1c = dpsi; C1c = dpsi-dpsi; D1c = 2*dpsi*dpsi;
            A2c = dpsi; B2c = dpsi; C2c = -(dpsi+dpsi); D2c = dpsi*dpsi^2+dpsi^2*dpsi;

            %VELOCITY__________________________________________
            dmurhou_dpsi(j,i) =( A1c*mu(j,i+1)*rho(j,i+1)*u(j,i+1) - B1c*mu(j,i-1)*rho(j,i-1)*u(j,i-1) -  C1c*mu(j,i)*rho(j,i)*u(j,i)) / (D1c);
            du_dpsi(j,i) = ( (A1c*(u(j,i+1))) - (B1c*(u(j,i-1))) - (C1c*(u(j,i))) )/(D1c );
            d2u_dpsi2(j,i) = 2*( (A2c*(u(j,i-1))) + (B2c*(u(j,i+1))) + (C2c*(u(j,i))) )/(D2c );
            
            %TEMPERATURE
            dkrhou_dpsi(j,i) = ( (A1c*k(j,i+1)*rho(j,i+1)*u(j,i+1)) - (B1c*(k(j,i-1)*rho(j,i-1)*u(j,i-1))) -  (C1c*(k(j,i)*rho(j,i)*u(j,i))))/(D1c);
            dT_dpsi(j,i) = ( (A1c*(T(j,i+1))) - (B1c*(T(j,i-1))) - (C1c*(T(j,i))) )/(D1c );
            d2T_dpsi2(j,i) = 2*( (A2c*(T(j,i-1))) + (B2c*(T(j,i+1))) + (C2c*(T(j,i))) )/(D2c);

         
        elseif i == gridlengthpsi(end)%Last_________________________________
            
           A1l = -1;           B1l = 0;        C1l = -1;        D1l = -dpsi;
           A2l = dpsi+dpsi;    B2l = -dpsi;    C2l = -dpsi;     D2l = (dpsi+dpsi)*dpsi^2-(dpsi+dpsi)^2*dpsi;

            %VELOCITY
            dmurhou_dpsi(j,i) = ( A1l*mu(j,i-1)*rho(j,i-1)*u(j,i-1)- 0 - C1l*(mu(j,i)*rho(j,i)*u(j,i)))/(D1l);
            du_dpsi(j,i) = ( A1l*u(j,i-1) - 0 - C1l*u(j,i))/(D1l);
            d2u_dpsi2(j,i) = 2*(A2l*u(j,i-1) + B2l*u(j,i-2) + C2l*u(j,i) )/(D2l)   ;
            %TEMPERATURE
             dkrhou_dpsi(j,i) = ( A1l*k(j,i-1)*rho(j,i-1)*u(j,i-1) - 0 - C1l*(k(j,i)*rho(j,i)*u(j,i)))/(D1l);
            dT_dpsi(j,i) =  ( A1l*T(j,i-1) - 0 - C1l*T(j,i))/(D1l);
            d2T_dpsi2(j,i) =  2*(A2l*T(j,i-1) + B2l*T(j,i-2) + C2l*T(j,i) )/(D2l)   ;
        end     
        
        
        if i == 1
            du = 0;
            dT = 0;
        
        elseif i>1
        % Chnage in velcoity
        du = (-1/(rho(j,i)*u(j,i))*dp_dx(j)+(dmurhou_dpsi(j,i))*du_dpsi(j,i)+mu(j,i)*rho(j,i)*u(j,i)*d2u_dpsi2(j,i))*dx(j);
        
        % Change in temperature
        dT = ((1/(rho(j,i)*cp(j,i)) * dp_dx(j) + mu(j,i)*rho(j,i)*u(j,i)/cp(j,i) * (du_dpsi(j,i))^2 + ...
            1/cp(j,i) * (dkrhou_dpsi(j,i) *dT_dpsi(j,i) + k(j,i)*rho(j,i)*u(j,i) * d2T_dpsi2(j,i)))*dx(j));
        end
        
        u(j+1,i) =u(j,i)+du;
        T(j+1,i) = T(j,i) + dT;
        
        % Dynamic viscosity
        mu(j+1,i) = sqrt(T(j+1,i)).*((1.1E-21).*T(j+1,i).^5 - (4.7E-18).*T(j+1,i).^4+(7.6E-15).*T(j+1,i).^3-(6.3E-12).*T(j+1,i).^2+(3E-9).*T(j+1,i)+5.6E-7);
        %Initialize
        rho(j+1,i) = P(j)./(R.*T(j+1,i));
        % Cp for air function
        cp(j+1,i) = (7.9e-17*T(j+1,i).^5 - 1.2e-013*T(j+1,i).^4 - 3.6e-10*T(j+1,i).^3 + 8.4e-7*T(j+1,i).^2 - 0.00035.*T(j+1,i) + 1)*1000;
        % K for air function
        k(j+1,i) = 1.7e-11*T(j+1,i).^3 - 5.1e-8*T(j+1,i).^2 + 0.0001*T(j+1,i) - 0.00053;
        %Mach number
        a(j+1,i) = sqrt(abs(gamma * R * T(j+1,i)));
        M(j+1,i) = u(j+1,i)./a(j+1,i);
   
    end

%     figure(20)
%     plot(T(j+1,:),'linewidth',2)
%     hold on
%   %  ylim([450 2100])
%     figure(11)
%      plot(u(j+1,:),'linewidth',2)
%      hold on

    if abs(x(j)-.1/100)< 1e-4
        location01  = j;
    elseif abs(x(j)-.2/100)< 1e-4
        location02  = j;
    end


    ydatathis(j,:) = (psi./(u(j+1,:).*rho(j+1,:)));
    nozzshape(j) = (psi(end)./(u(j+1,end).*rho(j+1,end))).*100./2.54;
end

ydatathis(:,1) = 0;

%% BOUNDARY LAYER PROFILES


%VELOCITY
uplothere = u;
%INLET
figure
subplot(1,5,1)
plot(uplothere(1,:),eta.*delta,'linewidth',2); grid on
hXLabel=xlabel('$u$ (m/s)','interpreter','latex'); hYLabel=ylabel('y (cm) ','interpreter','latex');
set(gca, 'FontName', 'Times'); set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'LineWidth', 1, 'fontsize',12)
htitle=title('$\textbf{Inlet (x = 0 \ cm)}$','interpreter','latex');
%THROAT
subplot(1,5,2)
plot(uplothere(338,:),ydatathis(338,:)*100,'linewidth',2); grid on
hXLabel=xlabel('$u$ (m/s)','interpreter','latex'); hYLabel=ylabel('y (cm) ','interpreter','latex');
set(gca, 'FontName', 'Times'); set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'LineWidth', 1, 'fontsize',12)
htitle=title('$\textbf{Throat (x = 0.037)}$','interpreter','latex');
%X = 0.1 cm
subplot(1,5,3)
plot(uplothere(location01,:),ydatathis(location01,:)*100,'linewidth',2); grid on
hXLabel=xlabel('$u$ (m/s)','interpreter','latex'); hYLabel=ylabel('y (cm) ','interpreter','latex');
set(gca, 'FontName', 'Times'); set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'LineWidth', 1, 'fontsize',12)
htitle=title('$\textbf{(x = 0.1 \ cm)}$','interpreter','latex');
%X = 0.2 cm
subplot(1,5,4)
plot(uplothere(location02,:),ydatathis(location02,:)*100,'linewidth',2); grid on
hXLabel=xlabel('$u$ (m/s)','interpreter','latex'); hYLabel=ylabel('y (cm) ','interpreter','latex');
set(gca, 'FontName', 'Times'); set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'LineWidth', 1, 'fontsize',12)
htitle=title('$\textbf{(x = 0.2 \ cm)}$','interpreter','latex');
%Exit
subplot(1,5,5)
plot(uplothere(gridlengthx,:),ydatathis(gridlengthx,:)*100,'linewidth',2); grid on
hXLabel=xlabel('$u$ (m/s)','interpreter','latex'); hYLabel=ylabel('y (cm) ','interpreter','latex');
set(gca, 'FontName', 'Times'); set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'LineWidth', 1, 'fontsize',12)
htitle=title('$\textbf{Exit (x = 0.3 \ cm)}$','interpreter','latex');



%TEMPERATURE
tplothere = T;
%INLET
figure;
subplot(1,5,1)
plot(tplothere(1,:),eta.*delta,'linewidth',2); grid on
hXLabel=xlabel('$T$ (K)','interpreter','latex'); hYLabel=ylabel('y (cm) ','interpreter','latex');
set(gca, 'FontName', 'Times'); set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'LineWidth', 1, 'fontsize',12)
htitle=title('$\textbf{Inlet (x = 0 \ cm)}$','interpreter','latex');
%THROAT
subplot(1,5,2)
plot(tplothere(338,:),ydatathis(338,:),'linewidth',2); grid on
hXLabel=xlabel('$T$ (K)','interpreter','latex'); hYLabel=ylabel('y (cm) ','interpreter','latex');
set(gca, 'FontName', 'Times'); set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'LineWidth', 1, 'fontsize',12)
htitle=title('$\textbf{Throat (x = 0.037)}$','interpreter','latex');
%X = 0.1 cm
subplot(1,5,3)
plot(tplothere(location01,:),ydatathis(location01,:),'linewidth',2); grid on
hXLabel=xlabel('$T$ (K)','interpreter','latex'); hYLabel=ylabel('y (cm) ','interpreter','latex');
set(gca, 'FontName', 'Times'); set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'LineWidth', 1, 'fontsize',12)
htitle=title('$\textbf{(x = 0.1 \ cm)}$','interpreter','latex');
%X = 0.2 cm
subplot(1,5,4)
plot(tplothere(location02,:),ydatathis(location02,:),'linewidth',2); grid on
hXLabel=xlabel('$T$ (K)','interpreter','latex'); hYLabel=ylabel('y (cm) ','interpreter','latex');
set(gca, 'FontName', 'Times'); set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'LineWidth', 1, 'fontsize',12)
htitle=title('$\textbf{(x = 0.2 \ cm)}$','interpreter','latex');
%EXIT
subplot(1,5,5)
plot(tplothere(gridlengthx,:),ydatathis(gridlengthx,:),'linewidth',2); grid on
hXLabel=xlabel('$T$ (K)','interpreter','latex'); hYLabel=ylabel('y (cm) ','interpreter','latex');
set(gca, 'FontName', 'Times'); set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'LineWidth', 1, 'fontsize',12)
htitle=title('$\textbf{Exit (x = 0.3 \ cm)}$','interpreter','latex');


%MACH


%% NOZZLE SHAPE
[x1,y1] = size(ydatathis);

xdatathis = ones(x1,y1);
for i = 1:y1
xdatathis(:,i) = x.*100;
end

%Fix the data matrices
xdatathis = xdatathis';
ydatathis = ydatathis';
u = flip(u(2:end,:)');
M = flip(M(2:end,:)');
T = flip(T(2:end,:)');




%Pressure Profile
figure
xplot2 = smooth(smooth(smooth(smooth(xplot)))*100)'; P_Po2 = smooth(smooth(smooth(smooth(P_Po2))))';
plot(smooth(smooth(xplot2)),smooth(smooth(P_Po2)),'k-','linewidth',3)
hold on
plot([-.4 0],[1 1],'k-','linewidth',3)
grid on
hXLabel=ylabel('$P/P_0$','interpreter','latex');
hYLabel=xlabel('Streamwise Distance (cm) ','interpreter','latex');
set(gca, 'FontName', 'Times')
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'LineWidth', 1, 'fontsize',20)



% Nozzle and Pressure
figure
subplot(2,1,1)
plot((smooth(x)*100),smooth(nozzshape)*.06/nozzshape(1),'k','linewidth',3)
hold on
plot(smooth(x)*100,-smooth(nozzshape)*.06/nozzshape(1),'k','linewidth',3)
plot([0, -.4],[smooth(nozzshape(1))*.06/nozzshape(1) smooth(nozzshape(1))*.06/nozzshape(1)],'k','linewidth',3)
plot([0, -.4],[-smooth(nozzshape(1))*.06/nozzshape(1) -smooth(nozzshape(1))*.06/nozzshape(1)],'k','linewidth',3)
hXLabel=xlabel('Streamwise Distance (cm)','interpreter','latex');
hYLabel=ylabel('Spanwise Distance(cm)','interpreter','latex');
set(gca, 'FontName', 'Times')
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'LineWidth', 1, 'fontsize',20)
grid on
grid minor
axis equal
xlim([-.4 .3])

subplot(2,1,2)
xplot2 = smooth(smooth(smooth(smooth(xplot)))*100)'; P_Po2 = smooth(smooth(smooth(smooth(P_Po2))))';
plot(smooth(smooth(xplot2)),smooth(smooth(P_Po2)),'k-','linewidth',3)
hold on
plot([-.4 0],[1 1],'k-','linewidth',3)
grid on
hXLabel=ylabel('$P/P_0$','interpreter','latex');
hYLabel=xlabel('Streamwise Distance (cm) ','interpreter','latex');
set(gca, 'FontName', 'Times')
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'LineWidth', 1, 'fontsize',20)

%Nozzle contours and psi constant grid spacing
figure
subplot(2,1,2)
for i = 1:5:(51)
plot(xdatathis,ydatathis(i,:)*100,'k','linewidth',2)
hold on
end
hXLabel=xlabel('x (cm)','interpreter','latex');
hYLabel=ylabel('y(cm)','interpreter','latex');
set(gca, 'FontName', 'Times')
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'LineWidth', 1, 'fontsize',20)
grid on
axis equal
xlim([0 .3])
ylim([0 .12])

subplot(2,1,1)
for i = 1:50:( 3065)
xline(x(i)*100,'k','linewidth',2)
hold on
end
for i = 1:2:(51)
yline(psi(i),'k','linewidth',2)
hold on
end
hXLabel=xlabel('x (cm)','interpreter','latex');
hYLabel=ylabel('$\psi$ ','interpreter','latex');
set(gca, 'FontName', 'Times')
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'LineWidth', 1, 'fontsize',20)
grid on
xlim([0 .3])


%% Contour Plots

%Mach contour level
Mlevel = [(.1:0.2:.9),1,(1.2:.2:2)];

%MACH
figure
[C,h] =contourf((xdatathis),(ydatathis)*100,M,10,'LevelList',Mlevel,'ShowText','on');
hold on
clabel(C,h,'FontSize',15,'Color','k')
contourf((xdatathis),(-ydatathis*100),M,10,'LevelList',Mlevel)
yline(0,'-.k','linewidth',2)
h=colorbar;
grid minor
colormap(jet);
grid on
axis equal
%________________________________________________
ax = gca; 
set(get(h,'label'),'string','Mach Number','Interpreter','latex' ,'FontName', 'Times', 'FontSize', 16);
hXLabel=xlabel('x (cm)');
hYLabel=ylabel('y (cm)');
%-------------------------------------------------
set(gca, 'FontName', 'Times')
set([ hYLabel], 'FontName', 'Times')
set([hXLabel, hYLabel], 'FontSize', 16)

set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'LineWidth', 1,'fontsize',20) 



% figure
% pcolor((xdatathis),(ydatathis),M)
% hold on
% pcolor((xdatathis),(-ydatathis),M)
% yline(0,'-.k','linewidth',2)
% colorbar
% shading interp
% colormap(jet);
% grid minor
% grid on
% x0=10;
% y0=10;
% width=700*1.25;
% height=400*1.25;
% set(gcf,'position',[x0,y0,width,height])
% set(gcf,'InvertHardcopy','off')


 
%TEMPERATURE
figure
contourf((xdatathis),(ydatathis*100),T,20)
hold on
%contourf((xdatathis),(-ydatathis*100),T,20)
yline(0,'-.k','linewidth',2)
h=colorbar;

grid minor
colormap(jet);
grid on
axis equal
%________________________________________________
ax = gca; 
set(get(h,'label'),'string','Temperature K','Interpreter','latex' ,'FontName', 'Times', 'FontSize', 16);
hXLabel=xlabel('x (cm)');
hYLabel=ylabel('y (cm)');
%-------------------------------------------------
set(gca, 'FontName', 'Times')
set([ hYLabel], 'FontName', 'Times')
set([hXLabel, hYLabel], 'FontSize', 16)

set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'LineWidth', 1,'fontsize',20) 

% figure
% pcolor((xdatathis),(ydatathis),T)
% hold on
% pcolor((xdatathis),(-ydatathis),T)
% yline(0,'-.k','linewidth',2)
% colorbar
% colormap(jet);
% shading interp
% colormap(jet);
% grid minor
% grid on
% x0=10;
% y0=10;
% width=700*1.25;
% height=400*1.25;
% set(gcf,'position',[x0,y0,width,height])
% set(gcf,'InvertHardcopy','off')


%VELOCITY
figure
contourf((xdatathis),(ydatathis*100),u,20)
hold on
%contourf((xdatathis),(-ydatathis*100),u,20)
yline(0,'-.k','linewidth',2)
h=colorbar;

grid minor
colormap(jet);
grid on

axis equal
%________________________________________________
ax = gca; 
set(get(h,'label'),'string','Velocity (m/s)','Interpreter','latex' ,'FontName', 'Times', 'FontSize', 16);
hXLabel=xlabel('x (cm)');
hYLabel=ylabel('y (cm)');
%-------------------------------------------------
set(gca, 'FontName', 'Times')
set([ hYLabel], 'FontName', 'Times')
set([hXLabel, hYLabel], 'FontSize', 16)

set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'LineWidth', 1,'fontsize',20) 
% figure
% pcolor((xdatathis),(ydatathis),u)
% hold on
% pcolor((xdatathis),(-ydatathis),u)
% yline(0,'-.k','linewidth',2)
% colorbar
% shading interp
% colormap(jet);
% grid minor
% grid on
% x0=10;
% y0=10;
% width=700*1.25;
% height=400*1.25;
% set(gcf,'position',[x0,y0,width,height])
% set(gcf,'InvertHardcopy','off')



%% NOW TUBE

%Create the psi space for this entire code
y_centerline = 0.0006/100
psi_centerline = y_centerline*V_in*rho_in;
dpsi = psi_centerline/50;
psi = (0:dpsi:psi_centerline);

xin = linspace(0.00001,0.3/100,300);
for j = 1:length(xin)
% INITIAL BOUNDARY LAYER___________________________________________________
% Compute Boundary layer thickness
delta(j) = 5.83*sqrt(nu_in*xin(j)/V_in);
%Relationship to blasius boundary layer
Capdelta = Pr_air^(1/3);
delta_T = delta(j)*Capdelta;

%Initialize Variables
y_from_psi(j,1) = 0;         u_Uinfty(1) = 0;        eta(1) = 0;   u_cube(1,j) = 0;
rho(1) = rho_wall;        TempNON_BL(1) = 1;      eta_T(1) = 0;
T_cube(j,1) = T_wall;            
for i = 2:length(psi)

    if i ==2
        y_from_psi(j,i) = y_from_psi(j,i-1) + dpsi/(rho(i-1)*51);
    else
        y_from_psi(j,i) = y_from_psi(j,i-1) + dpsi/(rho(i-1)*u_cube(i-1,j));
    end
    % BLASIUS BOUNDARY LAYER_______________________________________________
    % Nondim parameter to Blasius solution
    eta(i) = y_from_psi(j,i)./delta(j);
    % Blaius polynomial solution
    u_Uinfty(i) = 2*eta(i) - 2* eta(i).^3 + eta(i).^4;
    u_cube(i,j) = u_Uinfty(i)*V_in;
    %TEMPERATURE BOUNDARY PROFILE__________________________________________
    % Nondim parameter to temperaure solution
    eta_T(i) = y_from_psi(j,i)./delta_T;
    % Nondimesnional temperature expression
    TempNON_BL = 1 - 2*eta_T(i) + 2 * eta_T(i).^3 - eta_T(i).^4;
    T_cube(i,j) = TempNON_BL*(T_wall - T_in)+T_in;
    if u_cube(i,j) > V_in
        u_cube(i,j) = V_in;
    end
    if T_cube(i,j) > T_in
        T_cube(i,j) = T_in;
    end
    %Desnity change
    rho(i) = P_in/(R*T(i));
    a_cube = sqrt(gamma *R*T_cube(i,j));
    M_cube(i,j) = u_cube(i,j)./a_cube;
end

end

[meshx,meshy] = meshgrid(-xin*100,y_from_psi(1,:)*100);

%Tube SOLO
% figure
% contourf((meshx),(meshy),flip(u_cube),10);
% colorbar

xtotal = [meshx, xdatathis];
ytotal = [meshy, ydatathis];
utotal = [flip(u_cube),u];
Ttotal = [flip(T_cube(1:51,1:end)),T];
Mtotal = [flip(M_cube),M];

%% Contours with all

%VELOCITY
figure
[C,h] =contourf((xtotal),(ytotal),utotal,10);

hold on
contourf((xtotal),(-ytotal),utotal,10)
yline(0,'-.k','linewidth',2)
colorbar
grid minor
colormap(jet);
grid on
x0=10;
y0=10;
width=700*1.25;
height=400*1.25;
set(gcf,'position',[x0,y0,width,height])
set(gcf,'InvertHardcopy','off')


%MACH
figure
[C,h] =contourf((xtotal),(ytotal),Mtotal,10);
hold on
contourf((xtotal),(-ytotal),Mtotal,10)
yline(0,'-.k','linewidth',2)
colorbar
grid minor
colormap(jet);
grid on
x0=10;
y0=10;
width=700*1.25;
height=400*1.25;
set(gcf,'position',[x0,y0,width,height])
set(gcf,'InvertHardcopy','off')


% TEMPERATURE
figure
[C,h] =contourf((xtotal),(ytotal),Ttotal,10);
hold on
contourf((xtotal),(-ytotal),Ttotal,10)
yline(0,'-.k','linewidth',2)
colorbar
grid minor
colormap(jet);
grid on
x0=10;
y0=10;
width=700*1.25;
height=400*1.25;
set(gcf,'position',[x0,y0,width,height])
set(gcf,'InvertHardcopy','off')




