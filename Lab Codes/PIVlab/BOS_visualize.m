%% Load frames
Mconversion = 1;
col1 = 1; col2 = size(cell2mat(x(1)),2);
row1 = 1; row2 = size(cell2mat(x(1)),1);
clear u uL vL MAG rhoxL rhoyL MAG MAG_ALL NO_MEAN Current_uL
for FRAME = 1:size(v_original,1)
Displacement_Poisson.x = cell2mat(x(FRAME)).*Mconversion;
Displacement_Poisson.y = cell2mat(y(FRAME)).*Mconversion;
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

%% Create X-Displacement Movie

tic;
close all;
t_step = 0.05; % (s) time between frames
t_min = 0; t_max = 25;
vidname = 'X-Displacement_T#';
t_Fig = t_min:t_step:t_max;
writerObj = VideoWriter(['C:\Users\dlawson6\Desktop\Movies', vidname]);
writerObj.FrameRate = round(1/t_step);

open(writerObj);
figure('paperpositionmode', 'auto', 'position', [500 0 1024 1024],'DefaultAxesFontSize',30)
% Change above brackets accordingly: ['fig position from left of screen' 'fig position from bottomn of screen' 'Resolution X' 'Resolution Y']

t_int = (2200 / ((length(x) - 1))) / 1000;

for i = 1:length(x)
% i = 1:[Number of frames you have data for]
    i;

% All the normal figure stuff is below
colormap('jet')
c = colorbar;
f = 40;

c.Location = 'eastoutside';
c.Label.String = "X-Displacement";

xdisp = (flip(uL(:,:,i)));
pcolor(xdisp),shading interp,axis equal, axis tight

%set(gca,'color','k',FontSize=f);
xlabel('Pixel Index',FontSize=f);
ylabel('Pixel Index',FontSize=f);

xticks([20:20:160]);
yticks([20:20:160]);

%ylabel(c,'X-Displacement',FontSize=f);

% Set video background color

%set(gcf,'color','w');

% Insert Running Time Scale, WARNING: i does not appear to coincide with
% time like expected. Instead, use t_int = 2200 milliseconds / (length(x) - 1)

title(strcat('t = ', num2str(0.0+(i-1)*t_int),' s'), 'fontsize', 25)

% Write Video file

writeVideo(writerObj, getframe(gcf));

end
close(writerObj);
toc;

%% Plot Displacement in X, Y, or Magnitude

close all

frame = 700;
MAG = (flip(MAG_ALL(:,:,frame)));
xdisp = (flip(uL(:,:,frame)));
ydisp = (flip(vL(:,:,frame)));

f = 40;

%xdisp
figure
pcolor(xdisp),shading interp,axis equal, axis tight
c = colorbar;
set(gca,'color','k',FontSize=f)
%title('X-Displacement',FontSize=24);
ylabel('Pixel Index',FontSize=f);
xlabel('Pixel Index',FontSize=f);
xticks([20:20:160]);
yticks([20:20:160]);

ylabel(c,'X-Displacement',FontSize=f)


%ydisp
figure
pcolor(ydisp),shading interp,axis equal, axis tight
c = colorbar;
set(gca,'color','k',FontSize=f)
%title('Y-Displacement',FontSize=f);
ylabel('Pixel Index',FontSize=f);
xlabel('Pixel Index',FontSize=f);
xticks([20:20:160]);
yticks([20:20:160]);

ylabel(c,'Y-Displacement', FontSize=24)


%MAG
figure
pcolor(MAG),shading interp,axis equal, axis tight
c = colorbar;
set(gca,'color','k',FontSize=f)
%title('Displacement Magnitude',FontSize=f);
ylabel('Pixel Index',FontSize=f);
xlabel('Pixel Index',FontSize=f);
xticks([20:20:160]);
yticks([20:20:160]);

ylabel(c,'Displacement Magnitude', FontSize=f)
