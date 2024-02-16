function Data=pnrf_converter_f(blocks)
 
FromDisk = actxserver('Perception.Loaders.pNRF');
for loopd=1:length(blocks)
    clear Perception_Raw
    name=sprintf(blocks{loopd});
    theFiles = dir([name '*.pNRF']);
    for loopb=1:length(theFiles)
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
        save([blocks{loopd} '.mat'],'Perception_Raw','-v7.3')
    else
        fprintf(['Variable Perception_Raw not created for' name '\n'])
    end
end
end