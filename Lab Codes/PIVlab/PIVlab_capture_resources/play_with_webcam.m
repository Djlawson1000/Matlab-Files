clear
cam = webcam(1);
disp('Using the first webcam from the list of available webcams.')
%% Setting camera to a suitable resolution
desired_resolution = 1280; %this is the desired resolution for this example
available_resolutions=cam.AvailableResolutions;
width=zeros(numel(available_resolutions));
for i=1:numel(available_resolutions)
    temp = strsplit(available_resolutions{i},'x');
    width(i)=str2double(temp{1});
end
index_of_suitable_resolution = find(width>desired_resolution-100 & width < desired_resolution+100);
if isempty(index_of_suitable_resolution)
    index_of_suitable_resolution = numel(available_resolutions); %take the last listed resolution
end
index_of_suitable_resolution=index_of_suitable_resolution(1); %take only the first element if more suitable resolutions were found
cam.Resolution=available_resolutions{index_of_suitable_resolution};

%% Setting other camera parameters (adjust if needed)
cam.Exposure=-6;
cam.Contrast=40;
cam.Sharpness=99;
%% Prepare display
A=snapshot(cam);
A=rgb2gray(A);

img=snapshot(cam);
ima=imagesc(img);axis image;

start_acquisition=tic;
i=0;
while i < 50
	img=snapshot(cam);
	ima.CData=img;
	drawnow limitrate
	i=i+1;
end
framerate=i/toc(start_acquisition);
