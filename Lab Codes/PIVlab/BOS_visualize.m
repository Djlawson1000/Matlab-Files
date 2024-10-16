%% Load frames
Mconversion = 2.74; % was 0.71 length/Pixel Factor, change per model and adjust tick marks to account for pixel count
col1 = 1; col2 = size(cell2mat(x(1)),2);
row1 = 1; row2 = size(cell2mat(x(1)),1);
clear u uL vL MAG rhoxL rhoyL MAG MAG_ALL NO_MEAN Current_uL
for FRAME = 1:size(v_original,1)
% Displacement_Poisson.x = cell2mat(x(FRAME)).*Mconversion;
% Displacement_Poisson.y = cell2mat(y(FRAME)).*Mconversion;
Displacement_Poisson.u = cell2mat(u_original(FRAME)).*Mconversion.*1;
Displacement_Poisson.v =cell2mat(v_original(FRAME)).*Mconversion.*1;
[Disp] = (Displacement_Poisson);
uL(:,:,FRAME) = Disp.u(row1:row2,col1:col2);
vL(:,:,FRAME) = Disp.v(row1:row2,col1:col2);
nan_indices = isnan(uL(:,:,FRAME));
Current_uL = uL(:,:,FRAME);
Current_vL = vL(:,:,FRAME);
MAG_ALL(:,:,FRAME) = sqrt(Current_uL.^2 + Current_vL.^2);
%Output
if mod(FRAME, 100) == 0
disp(['Iteration ', num2str(FRAME)]);
end
end
Displacement_Poisson.x = meshgrid(1:size(uL(:,:,FRAME),2),1:size(uL(:,:,FRAME),1))*Mconversion;
Displacement_Poisson.y =meshgrid(1:size(uL(:,:,FRAME),1),1:size(uL(:,:,FRAME),2))'*Mconversion;

%% Create X-Displacement Movie

tic;
close all;
t_step = 0.05; % (s) time between frames
t_min = 0; t_max = 25;
vidname = 'X-Displacement_T35';
t_Fig = t_min:t_step:t_max;
writerObj = VideoWriter(['C:\Users\dlawson6\Desktop\Movies', vidname]);
writerObj.FrameRate = round(1/t_step);

open(writerObj);
figure('paperpositionmode', 'auto', 'position', [500 0 1200 1200],'DefaultAxesFontSize',30)
% Change above brackets accordingly: ['fig position from left of screen' 'fig position from bottomn of screen' 'Resolution X' 'Resolution Y']

t_int = (2200 / ((length(x) - 1))) / 1000;

for i = 1:length(x)
% i = 1:[Number of frames you have data for]
    i;

% All the normal figure stuff is below
f = 40;

xdisp = (flip(uL(:,:,i)));
pcolor(Displacement_Poisson.x, Displacement_Poisson.y, xdisp),shading interp,axis equal, axis tight

xlabel('X (mm)',FontSize=f);
ylabel('Y (mm)',FontSize=f);

% xticks([20:20:160]);
% yticks([20:20:160]);

% Set video background color

set(gcf,'color','w');
set(gca,'color','k');


% Insert Running Time Scale, WARNING: i does not appear to coincide with
% time like expected. Instead, use t_int = 2200 milliseconds / (length(x) - 1)

title(strcat('t = ', num2str(0.0+(i-1)*t_int),' s'), 'fontsize', 25)

colorbar
colormap('hot')
c = colorbar;
ylabel(c,'X-Displacement (mm)',FontSize=f);
clim([-8 3]) % Adjust this based on the range seen during steady flow of each test

% Write Video file

writeVideo(writerObj, getframe(gcf));

end
close(writerObj);
toc;
%% Create Y-Displacement Movie

tic;
close all;
t_step = 0.05; % (s) time between frames
t_min = 0; t_max = 25;
vidname = 'Y-Displacement_T35';
t_Fig = t_min:t_step:t_max;
writerObj = VideoWriter(['C:\Users\dlawson6\Desktop\Movies', vidname]);
writerObj.FrameRate = round(1/t_step);

open(writerObj);
figure('paperpositionmode', 'auto', 'position', [500 0 1200 1200],'DefaultAxesFontSize',30)
% Change above brackets accordingly: ['fig position from left of screen' 'fig position from bottom of screen' 'Resolution X' 'Resolution Y']

t_int = (2200 / ((length(y) - 1))) / 1000;

for i = 1:length(y)
% i = 1:[Number of frames you have data for]
    i;

% All the normal figure stuff is below
f = 40;

ydisp = (flip(vL(:,:,i)));
pcolor(Displacement_Poisson.x, Displacement_Poisson.y, ydisp),shading interp,axis equal, axis tight

xlabel('X (mm)',FontSize=f);
ylabel('Y (mm)',FontSize=f);

% xticks([20:20:160]);
% yticks([20:20:160]);

% Set video background color

set(gcf,'color','w');
set(gca,'color','k');


% Insert Running Time Scale, WARNING: i does not appear to coincide with
% time like expected. Instead, use t_int = 2200 milliseconds / (length(x) - 1)

title(strcat('t = ', num2str(0.0+(i-1)*t_int),' s'), 'fontsize', 25)

colorbar
colormap('hot')
c = colorbar;
ylabel(c,'Y-Displacement (mm)',FontSize=f);
clim([0 0.2]) % Adjust this based on the range seen during steady flow of each test

% Write Video file

writeVideo(writerObj, getframe(gcf));

end
close(writerObj);
toc;

%% Create Magnitude Movie

tic;
close all;
t_step = 0.05; % (s) time between frames
t_min = 0; t_max = 25;
vidname = 'Disp Magnitude_T35';
t_Fig = t_min:t_step:t_max;
writerObj = VideoWriter(['C:\Users\dlawson6\Desktop\Movies', vidname]);
writerObj.FrameRate = round(1/t_step);

open(writerObj);
figure('paperpositionmode', 'auto', 'position', [500 0 1200 1200],'DefaultAxesFontSize',30)
% Change above brackets accordingly: ['fig position from left of screen' 'fig position from bottomn of screen' 'Resolution X' 'Resolution Y']

t_int = (2200 / ((length(x) - 1))) / 1000;

for i = 1:length(x)
% i = 1:[Number of frames you have data for]
    i;

% All the normal figure stuff is below
f = 40;

MAG = (flip(MAG_ALL(:,:,i)));
pcolor(Displacement_Poisson.x, Displacement_Poisson.y, MAG),shading interp,axis equal, axis tight

xlabel('X (mm)',FontSize=f);
ylabel('Y (mm)',FontSize=f);

% xticks([20:20:160]);
% yticks([20:20:160]);

% Set video background color

set(gcf,'color','w');
set(gca,'color','k');


% Insert Running Time Scale, WARNING: i does not appear to coincide with
% time like expected. Instead, use t_int = 2200 milliseconds / (length(x) - 1)

title(strcat('t = ', num2str(0.0+(i-1)*t_int),' s'), 'fontsize', 25)

colorbar
colormap('gray')
c = colorbar;
ylabel(c,'Displacement Magnitude (mm)',FontSize=f);
clim([0 0.2]) % Adjust this based on the range seen during steady flow of each test

% Write Video file

writeVideo(writerObj, getframe(gcf));

end
close(writerObj);
toc;


%% Plot Displacement in X, Y, or Magnitude

close all

frame = 1250;
MAG = (flip(MAG_ALL(:,:,frame)));
xdisp = (flip(uL(:,:,frame)));
ydisp = (flip(vL(:,:,frame)));

f = 40;

%xdisp
figure
pcolor(Displacement_Poisson.x, Displacement_Poisson.y, xdisp),shading interp,axis equal, axis tight
c = colorbar;
%clim([0.1 0.5]);
set(gca,'color','k',FontSize=f)
set(gcf,'color','w');
%title('X-Displacement',FontSize=24);
xlabel('X (mm)',FontSize=f);
ylabel('Y (mm)',FontSize=f);
% xticks([20:20:160]);
% yticks([20:20:160]);

ylabel(c,'X-Displacement (mm)',FontSize=f)


%ydisp
figure
pcolor(Displacement_Poisson.x, Displacement_Poisson.y, ydisp),shading interp,axis equal, axis tight
c = colorbar;
%clim([0.1 0.5]);
set(gca,'color','k',FontSize=f)
set(gcf,'color','w');
%title('Y-Displacement',FontSize=f);
xlabel('X (mm)',FontSize=f);
ylabel('Y (mm)',FontSize=f);
% xticks([20:20:160]);
% yticks([20:20:160]);

ylabel(c,'Y-Displacement (mm)', FontSize=f)


%MAG
figure
pcolor(Displacement_Poisson.x, Displacement_Poisson.y, MAG),shading interp,axis equal, axis tight
c = colorbar;
%clim([0.1 0.5]);
set(gca,'color','k',FontSize=f)
set(gcf,'color','w');
%title('Displacement Magnitude',FontSize=f);
xlabel('X (mm)',FontSize=f);
ylabel('Y (mm)',FontSize=f);
% xticks([20:20:160]);
% yticks([20:20:160]);

ylabel(c,'Displacement Magnitude (mm)', FontSize=f)

%%    mm / Pixel for each Test Run

% T1:  1.02
% T2:  1.02
% T3:  2.00 
% T13: 0.71
% T14: -
% T15: -
% T17: -
% T18: 0.71
% T19: 0.763






















