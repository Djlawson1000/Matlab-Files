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
%% Plot
close all

frame = 700;
MAG = (flip(MAG_ALL(:,:,frame)));
xdisp = (flip(uL(:,:,frame)));
ydisp = (flip(vL(:,:,frame)));


%xdisp
figure
pcolor(xdisp),shading interp,axis equal, axis tight
colorbar
set(gca,'color','k')

%ydisp
figure
pcolor(ydisp),shading interp,axis equal, axis tight
colorbar
set(gca,'color','k')

%MAG
figure
pcolor(MAG),shading interp,axis equal, axis tight
colorbar
set(gca,'color','k')