clear all
clc
close all
imaqreset

%% Prepare camera
delete(imaqfind); %clears all previous videoinputs
try
    hwinf = imaqhwinfo;
    %imaqreset
catch
    disp('Error: Image Acquisition Toolbox not available!')
end
info = imaqhwinfo(hwinf.InstalledAdaptors{1});
if strcmp(info.AdaptorName,'gentl')
    disp('gentl adaptor found.')
else
    disp('ERROR: gentl adaptor not found. Please install the GenICam / GenTL support package from here:')
    disp('https://de.mathworks.com/matlabcentral/fileexchange/45180')
end
OPTRONIS_name = info.DeviceInfo.DeviceName;


OPTRONIS_supported_formats = info.DeviceInfo.SupportedFormats;

% select bitmode (some support 8, 10, 12 bits)
bitmode =10; %even if the camera accepts 10 bitmode, the returned data still is  8bit in matlab
OPTRONIS_vid = videoinput(info.AdaptorName,info.DeviceInfo.DeviceID,['Mono' sprintf('%0.0d',bitmode)]);
info.DeviceInfo.DefaultFormat = ['Mono' sprintf('%0.0d',bitmode)];
OPTRONIS_settings = get(OPTRONIS_vid);
OPTRONIS_settings.PreviewFullBitDepth='On';
OPTRONIS_vid.PreviewFullBitDepth='On';

%{
%% recording to disk
filelocation = "C:\Users\william\Desktop\PIV_data_test";
filename = "recording1.avi";
fullFilename = fullfile(filelocation, filename);

% Create and configure the video writer
logfile = VideoWriter(fullFilename, "Grayscale AVI");

% Configure the device to log to disk using the video writer
OPTRONIS_vid.LoggingMode = "disk";
OPTRONIS_vid.DiskLogger = logfile;
%}


%% set camera parameters for triggered acquisition

%%{
% EXTERNER TRIGGER 
%[status,cmdout] = system('cd C:\Program Files\Euresys\eGrabber\bin\x86_64\ && gentl script singleframe.js', '-echo'); %das muss ausgeführt werden bevor kamera angesprochen wird.

triggerconfig(OPTRONIS_vid, 'hardware','DeviceSpecific','DeviceSpecific');
OPTRONIS_settings.TriggerSource = 'SingleFrame';
%OPTRONIS_settings.Source.ExposureMode = 'Timed';
OPTRONIS_settings.TriggerMode ='On';
%%}

%{
%INTERNER TRIGGER 
[status,cmdout] = system('cd C:\Program Files\Euresys\eGrabber\bin\x86_64\ && gentl script continuous.js', '-echo');
triggerconfig(OPTRONIS_vid, 'immediate','none','none');
OPTRONIS_settings.TriggerSource = 'Continuous';
%OPTRONIS_settings.Source.ExposureMode = 'Timed';
OPTRONIS_settings.TriggerMode ='Off';
%}

OPTRONIS_src=getselectedsource(OPTRONIS_vid);



%set(OPTRONIS_src, "AcquisitionMode", "SingleFrame"); get nicht

framerate=1700



if contains(OPTRONIS_name,'Cyclone-2-2000-M')
	disp(['Found camera: ' 'Cyclone-2-2000-M'])
    %framerate=floor(1/((expotime+3)/1000^2)) %muss man auch setzen damit exposure time akzeptiert wird...
    expotime=ceil(1/framerate*1000^2-3) %minimum ist 459 bei 2166 fps
    minexpo=2;
    warning('off','imaq:gentl:hardwareTriggerTriggerModeOff'); %trigger property of OPTRONIS cannot be set in Matlab.
    warning('off','MATLAB:JavaEDTAutoDelegation'); %strange warning
elseif contains (OPTRONIS_name,'Cyclone-1HS-3500-M')
	disp(['Found camera: ' 'Cyclone-1HS-3500-M'])
	%framerate=floor(1/((expotime+3)/1000^2)) %muss man auch setzen damit exposure time akzeptiert wird...
    expotime=ceil(1/framerate*1000^2-3) %3178 ist maximum
    minexpo=2;
    warning('off','imaq:gentl:hardwareTriggerTriggerModeOff'); %trigger property of OPTRONIS cannot be set in Matlab.
    warning('off','MATLAB:JavaEDTAutoDelegation'); %strange warning
elseif contains (OPTRONIS_name,'Cyclone-25-150-M')
	disp(['Found camera: ' 'Cyclone-25-150-M'])
	%framerate=floor(1/((expotime+3)/1000^2)) %muss man auch setzen damit exposure time akzeptiert wird...
    expotime=ceil(1/framerate*1000^2-24) %149 ist maximum
    minexpo=12;
    warning('off','imaq:gentl:hardwareTriggerTriggerModeOff'); %trigger property of OPTRONIS cannot be set in Matlab.
    warning('off','MATLAB:JavaEDTAutoDelegation'); %strange warning
else
    disp('camera type unknown!')
end

OPTRONIS_src.ExposureTime = minexpo;
OPTRONIS_src.AcquisitionFrameRate = framerate;

OPTRONIS_vid.FramesPerTrigger = 1165;%1000; %scheinbar das max. bei 64 gb = 16000
flushdata(OPTRONIS_vid);


%OPTRONIS_vid.ROIPosition = [96,100,512,500]
%OPTRONIS_vid.ROIPosition = [0 0 1920 1080]
start(OPTRONIS_vid);

OPTRONIS_src.AcquisitionFrameRate = framerate; %muss man auch setzen damit exposure time akzeptiert wird...
OPTRONIS_src.ExposureTime=expotime;

%Nur wenn man nach dem start expotime setzt wird es respektiert.
prev=preview(OPTRONIS_vid);


[userview,systemview] = memory;
%mem_base_use = systemview.PhysicalMemory.Total/1024^3 - systemview.PhysicalMemory.Available/1024^3 %memory in GB without images in RAM
mem_free=systemview.PhysicalMemory.Available/1024^3 % vor capture abfragen --> bestimmt anzahl bilder.
current_ROI = OPTRONIS_vid.ROIPosition;
%mem_required = current_ROI(3)*current_ROI(4)*bitmode/8/1024^3 * OPTRONIS_vid.FramesPerTrigger*2

ram_reserve=1; %how much RAM (in GB) is needed for Matlab operations not including image capture
max_possible_images=floor((mem_free-ram_reserve) / (current_ROI(3)*current_ROI(4)*bitmode/8/1024^3 *2))

%vor dem start ausrechnen ob es reicht. es wird 2x so viel benötigt wie man bilder aufnimmt. getdata beraucht genau nochmal so viel.

while OPTRONIS_vid.FramesAcquired < (OPTRONIS_vid.FramesPerTrigger)
    pause (0.1)

	%disp(num2str(bla.MaxPossibleArrayBytes/2/2073600))

	%ein Bild ist 2073600 bytes

	
end
stop(OPTRONIS_vid);
stoppreview(OPTRONIS_vid)
closepreview(OPTRONIS_vid)



OPTRONIS_data = getdata(OPTRONIS_vid,OPTRONIS_vid.FramesPerTrigger);


%{
cntr=0;
bitmultiplicator=1;
ImagePath = 'C:\Users\william\Documents\MATLAB\datatest\';
tic
for image_save_number=1:2:size(OPTRONIS_data,4)
    imgA_path=fullfile(ImagePath,['PIVlab_' sprintf('%4.4d',cntr) '_A.tif']);
    imgB_path=fullfile(ImagePath,['PIVlab_' sprintf('%4.4d',cntr) '_B.tif']);
    imwrite(OPTRONIS_data(:,:,:,image_save_number)*bitmultiplicator,imgA_path,'compression','none'); %tif file saving seems to be the fastest method for saving data...
    imwrite(OPTRONIS_data(:,:,:,image_save_number+1)*bitmultiplicator,imgB_path,'compression','none');
    cntr=cntr+1;
    disp(['saving frame ' num2str(cntr)])
end
toc

%}



