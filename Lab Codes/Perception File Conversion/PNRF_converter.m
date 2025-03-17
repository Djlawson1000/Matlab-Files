%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PNRF_converter.m
% Save this code in a folder with pNRF files. 
% Set the 'blocks' cell to the name of the pNRF files you are looking to convert. 
% The code will find the counts of the files in numerical order. 
% It will also determine which cartridges are being used, in the order they are used and turned on
% during the experiment so double check which channels are turned on in Perception. 
% Finally, it will figure out which channels are used. It should sort that by number. 
% A lot of the core code for extracting the data from the pNRF file is
% modified from the following link. 
% You need to download a file from this link:
% https://www.hbm.com/en/7557/hardware-and-software-interfacing-with-genesis-high-speed/
% Somewhere in the middle of the page is the pNRF reader toolkit. 
% Read the manual and example code. 
% Install the toolkit, the last link on the above webpage.
% Then run this code.
% Author: Erik Hoberg <ehoberg@nd.edu>
% Version: 1.0
% Version Date: 03/03/2023
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUTS
%blocks: input the handle of your calibration files and test files here.
%Make sure 'handle001.PNRF' exists! or adjust line 47 to an existant file
%with that handle.
%
%OUTPUTS
% A .mat file saved in the same folder of the same name as the handle in blocks. Each file contains
% all the files counted for that handle in sequential order from 1--999. 
clc; clear; close all

blocks={'N'};%,'kcal','dcal','aftertest'
pnrf_converter(blocks);
 
function []=pnrf_converter(blocks)
 1
FromDisk = actxserver('Perception.Loaders.pNRF');
for loopd=1:length(blocks)
    loopd
    clear Perception_Raw
    name=sprintf(blocks{loopd});
    theFiles = dir([name '*.pNRF']);
    for loopb=1:length(theFiles)
        loopb
        RecordingName =[theFiles(loopb).name];
        if isfile(RecordingName)==1
            Data = FromDisk.LoadRecording(RecordingName);
            for counta=1:20 %Finding which cartridges to read
                TFa=isempty(Data.Recorders.Item(counta));
                if TFa==0
                    Recorder = Data.Recorders.Item(counta);%
                    Channels = Recorder.Channels;
                    for loopa=1:20 %Finding which channels used on a cartridge and saving them
                         TFbb=isempty(Channels.Item(loopa));
                         if TFbb==0
                         Channel = Channels.Item(loopa);
                            ItfData = Channel.DataSource(3);
                            Sweeps = ItfData.Sweeps;
                            dStartTime = Sweeps.StartTime;
                            dEndTime = Sweeps.EndTime;
                            SegmentsOfData = ItfData.Data(dStartTime, dEndTime);
                            Segment = SegmentsOfData.Item(1);
                            NumberOfSamples = Segment.NumberOfSamples;
                            WaveformData = Segment.Waveform(5, 1, NumberOfSamples, 1)';
                            TFb=any(WaveformData);
                         else
                             TFb=0;
                         end
                        if TFb==1
                            tEnd = Segment.StartTime+(NumberOfSamples - 1)*Segment.SampleInterval;
                            t = Segment.StartTime: Segment.SampleInterval:tEnd;
                            Perception_Raw(loopb,counta).time=t;
                            Perception_Raw(loopb,counta).channel(:,loopa)=WaveformData;
                        end
                    end
                end
            end
        else
            fprintf(['\n' RecordingName ' Not valid!\n'])
        end
    end
    if exist('Perception_Raw','var') == 1
        1
        save([blocks{loopd} '.mat'],'Perception_Raw','-v7.3')
    else
        2
        fprintf(['Variable Perception_Raw not created for' name '\n'])
    end
end
end
