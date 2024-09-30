%% Density filtered back projection


%% INPUTS
Test = 1;              % Test case #
dt   = .25;            % Theta increment deg
M = 1;                 % Magnification factor keep 1 (For testing)

Frame = 400;
v_x = uL(:,:,Frame);
v_y = vL(:,:,Frame);

% Defines theta for sinogram
theta = [0:dt:179];
N = length(theta);

% Grid size - for ploting and integration
AxyDENS = v_x;
%Flip to make axisymmetric image
% Total_NORMx = [ flip(AxyDENS(1:end,100:end-100)); AxyDENS(1:end,100:end-100)];
Total_NORMx =  (AxyDENS(:,:)); 
% Total_NORMx = [ flip(AxyDENS(:,:)); (AxyDENS(:,:))];
nan_indices = isnan(Total_NORMx(:,:));
figure
surfc(Total_NORMx,'EdgeColor','none');
view(2); axis equal

AxyDENS = v_y;
%Flip to make axisymmetric image
% Total_NORMy = [ flip(AxyDENS(1:end,100:end-100)); AxyDENS(1:end,100:end-100)];
Total_NORMy = (AxyDENS(:,:)); 
% Total_NORMy = [ flip(AxyDENS(:,:)); (AxyDENS(:,:))];
nan_indices = isnan(Total_NORMy(:,:));
figure
surfc(Total_NORMy,'EdgeColor','none');
view(2); axis equal

[I,J] = size(Total_NORMy);% Get number of elements - axi-symmetric data


%% compute x - filtered back projection
dndx = zeros(I,J); % Allocate memory for refraction index gradient

%nan_indices = isnan(Total_NORMx(:,:));
Total_NORMx(nan_indices) = 0;
%
f = waitbar(0,'1','Name','Computing dndx...');
t0 =  tic;
for j=1:J
    dp = Total_NORMx(:,j);    % Allocate data at station x=x(j);
    sg = repmat(dp,1,N);% Create sinogram
    % Compute inverse radon transform
%     idp = iradon(sg,theta, 'v5cubic', 'Shepp-Logan', 1, I); % Using Shepp-Logan filter 
    idp = iradon(sg,theta,I); % Using deafult filter
    dndx(:,j) = idp(:,round(I/2)); % allocate slice center line at x=x(j)
    % Execution time
    tend = toc(t0);
    tend = (J-j)/j*tend;
    waitbar(j/J,f,sprintf('Estimated ending time %s',datetime+seconds(tend))); % Waitbar
end
close(f); % close wait bar
dndx(nan_indices==1) = nan; % Set solid body region to nan

figure
surf(dndx,'edgecolor','none');view(2)
%caxis([-1*10^(-4) 1*10^(-4)])

dndxNN = flip(dndx(1:size(AxyDENS(:,:),1),1:size(AxyDENS(:,:),2)));
figure
surf(dndxNN,'edgecolor','none');view(2)
%caxis([-.5*10^(-4) .5*10^(-4)])

%% compute y-filtered backprojection
dndy = zeros(I,J); % Allocate memory for refraction index gradient

nan_indices = isnan(Total_NORMy(:,:));
Total_NORMy(nan_indices) = 0;
%
f = waitbar(0,'1','Name','Computing dndy...');
t0 =  tic;
for j=1:J
    dp = Total_NORMy(:,j);    % Allocate data at station x=x(j);
    sg = repmat(dp,1,N);% Create sinogram
    % Compute inverse radon transform
%     idp = iradon(sg,theta, 'v5cubic', 'Shepp-Logan', 1, I); % Using Shepp-Logan filter 
    idp = iradon(sg,theta,I); % Using deafult filter
    dndy(:,j) = idp(:,round(I/2)); % allocate slice center line at x=x(j)
    % Execution time
    ttotal = toc(t0);
    tend = (J-j)/j*ttotal;
    waitbar(j/J,f,sprintf('Estimated ending time %s',datetime+seconds(tend))); % Waitbar
end
close(f); % close wait bar
dndy(nan_indices==1) = nan; % Set solid body region to nan

figure
surf(dndy,'edgecolor','none');view(2)
%caxis([-1*10^(-4) 1*10^(-4)])

dndyNN = flip(dndy(1:size(AxyDENS(:,:),1),1:size(AxyDENS(:,:),2)));
figure
surf(dndyNN,'edgecolor','none');view(2)
%caxis([-1*10^(-4) 1*10^(-4)])

