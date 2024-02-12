clear all
clc
close all
imaqreset

imaqmex('feature','-debugGenTLDiscovery',true)
imaqmex('feature','-debugGenTLAcquisition',true)
imaqmex('feature','-debugGenTLOperation',true)

hwinf = imaqhwinfo;
info = imaqhwinfo(hwinf.InstalledAdaptors{1});

OPTRONIS_vid = videoinput(info.AdaptorName,info.DeviceInfo.DeviceID,'Mono10');

OPTRONIS_settings = get(OPTRONIS_vid);
OPTRONIS_settings.PreviewFullBitDepth='On';
OPTRONIS_vid.PreviewFullBitDepth='On';

OPTRONIS_src=getselectedsource(OPTRONIS_vid);

%Because it is not possible to configure the trigger in Matlab, I have set
%the camera to external triggering using the manufacturers software. I then
%supply a continuous hardware trigger signal to the camera.
%{ 
%not having any effect:
OPTRONIS_src.AcquisitionMode='SingleFrame';
triggerconfig(OPTRONIS_vid, 'hardware','DeviceSpecific','DeviceSpecific');
%OPTRONIS_settings.TriggerSource = 'SingleFrame';
OPTRONIS_settings.Source.ExposureMode = 'Timed';
%OPTRONIS_settings.TriggerMode ='On';
%}


framerate=20;


expotime=ceil(1/framerate*1000^2-3) %The camera supports e.g. 2166 fps with an exposure time of 459 µs max.
minexpo=2;

OPTRONIS_src.ExposureTime = minexpo; %exposure time needs to be set to the smallest possible value first, otherwise high framerates are not accepted
OPTRONIS_src.AcquisitionFrameRate = framerate;
OPTRONIS_vid.FramesPerTrigger = 1;



%OPTRONIS_vid.ROIPosition = [0 0 1920 1080];
OPTRONIS_vid.ROIPosition = [512,512,512,500]; %When setting this ROI, the y-offset will always be ignored. The camera reports the correct ROI in ROIPosition, but in reality, the ROI is always starting at the very upper edge of the sensor.

ROI=OPTRONIS_vid.ROIPosition %reports the correct numbers

start(OPTRONIS_vid);
%Only when frame rate and exposure time are set again AFTER starting the acquisition, exposure time is accepted. Otherwise it is ignored

OPTRONIS_src.AcquisitionFrameRate = framerate; 
OPTRONIS_src.ExposureTime=expotime;

prev=preview(OPTRONIS_vid);


while OPTRONIS_vid.FramesAcquired < (OPTRONIS_vid.FramesPerTrigger)
    pause (0.1)
end
stop(OPTRONIS_vid);
stoppreview(OPTRONIS_vid)
closepreview(OPTRONIS_vid)
OPTRONIS_data = getdata(OPTRONIS_vid,OPTRONIS_vid.FramesPerTrigger,'native'); %The returned variable will always be 4-D uint8, no matter if I set Mono8 or Mono10. Even if I force uint16, then the data is 8 bit in reality (can be determined by counting the amount of possible grayscale steps per pixel).

% disable debugging
pause(1)
imaqmex('feature','-debugGenTLDiscovery',false)
imaqmex('feature','-debugGenTLAcquisition',false)
imaqmex('feature','-debugGenTLOperation',false)