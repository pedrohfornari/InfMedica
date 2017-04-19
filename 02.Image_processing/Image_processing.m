%% Informatica Medica Trabralho 2
 % Processamento de Imagens
 % Pedro Henrique Fornari
 % 13104320

 clear all;
 close all;
%% Load Image
I = imread('bigben.png'); %load image
Igray = rgb2gray(I); %convert it to gray scale
figure ('Name', 'Original Image');
imshow(I);

%% Apply bit plane slicing
B = cell(8,1);

figure('Name', 'Bit Panel Slice')
for i = 1:8
    B{i}=zeros(size(Igray)); %Pre set each cell
    B{i}=bitset(B{i},i,bitget(Igray,i)); %get the selected bit
    B{i}=uint8(B{i}); %transform it to an 8 bit integer
    subplot(2, 4, i); 
    imshow(B{i}*255); %multiply by 255 to correct the preset
    str = sprintf('Panel Bit %d', i); %update title
    title(str);
end

%% Get RGB channels

Rchannel = I(:, :, 1); %Red channel is in the first dimension
Gchannel = I(:, :, 2); %Green channel is in the second dimension
Bchannel = I(:, :, 3); %Blue channel is in the third dimension

figure
subplot(1, 3, 1);
imshow(Rchannel);
title('Red Channel');
subplot(1, 3, 2);
imshow(Gchannel);
title('Green Channel');
subplot(1, 3, 3);
imshow(Bchannel);
title('Blue Channel');

%% Get gray image by manual conversion

% define each value
alfa = 0.2989; beta = 0.5870; gamma = 0.1140;

% transform each pixel to gray scale
for i = 1:size(I, 1)
    for j = 1:size(I, 2)
        Imgray(i, j) = alfa*I(i, j, 1) + beta*I(i, j, 2) + gamma*I(i, j, 3);
    end
end

figure('Name', 'Manual conversion from rgb to gray scale')
imshow(Imgray);

%% Get hue, saturation and value from a RGB image

HSVImage = rgb2hsv(I);

Hue = HSVImage(:,:,1);
Saturation = HSVImage(:,:,2);
Value = HSVImage(:,:,3);


figure
subplot(1, 3, 1);
imshow(Hue);
title('Hue');
subplot(1, 3, 2);
imshow(Saturation);
title('Saturation');
subplot(1, 3, 3);
imshow(Value);
title('Value');

%% Do some adition operations with images

%sum 100 to an image
Iplus100 = imadd(I, 100);

%subtract an image from another one
A = imread('bigben.png');
B = imread('edin_castle.png');
AminusB = imsubtract(A,B);

%invert image
inverted = imadd(-A, 255);    

%tresholding image
