cd C:\Users\dlawson6\Desktop\'Matlab Files'\'Lab Codes'

%cd /Users/dan/Documents/GitHub/'Matlab-Files'/'Lab Codes'

clear;

%Folder where you have the Tiff stacks
cd 'TiffStack Conversion'

% Name of Wind on Tiff stack
Data = tiffreadVolume('T70_1101_2201.tif');

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
        imwrite(I,['_110' int2str(i), '.tiff']);
    elseif i <100
        imwrite(I,['_11' int2str(i), '.tiff']);
    elseif i <1000
        imwrite(I,['_1' int2str(i), '.tiff']);
    elseif i< 10000
        imwrite(I,['_' int2str(i), '.tiff']);
    end

    disp(i)
end

%% REFERENCE

%Name of Wind off (Reference) Tiff Stack
Data = tiffreadVolume('T47Ref.tif');

i = 1;
I = Data(:,:,i);
if Crop ==1
    I = imcrop(I,rect);
end
imwrite(I,['0000_Ref','.tiff']);
