a = imread('78_0001.tiff');
b = imread('80_0001.tiff');

c = a - b;

%% 

imshow(c)