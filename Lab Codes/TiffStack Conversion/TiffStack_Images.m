cd /Users/jonathandavami/Documents/MATLAB/AAA_Hypersonic_Research/ROTEXT_Jan_2024/BOS/Saved_Frames/FramesPhantom/

clear;

%Folder where you have the Tiff stacks
cd T29Frames

% Name of Wind on Tiff stack
Data = tiffreadVolume('T1.tif');

%% WIND ON

%Crop images? Yes then set to 1 || Or No set to 0
Crop = 0;

numFrames= size(Data,3);
n=numFrames;
for i = 1:1:n
    I = Data(:,:,i);

    if Crop ==1
        if i == 1
            figure
            [J,rect] = imcrop(I);
            I = imcrop(I,rect);
        else
            I = imcrop(I,rect);
        end
    end


    if i <10
        imwrite(I,['_000' int2str(i), '.tiff']);
    elseif i <100
        imwrite(I,['_00' int2str(i), '.tiff']);
    elseif i <1000
        imwrite(I,['_0' int2str(i), '.tiff']);
    elseif i< 10000
        imwrite(I,['_' int2str(i), '.tiff']);
    end

    disp(i)
end

%% REFERENCE

%Name of Wind off (Reference) Tiff Stack
Data = tiffreadVolume('T29_1_2.tif');

i = 1;
I = Data(:,:,i);
if Crop ==1
    I = imcrop(I,rect);
end
imwrite(I,['0000_Ref','.tiff']);
