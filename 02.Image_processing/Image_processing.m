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
figure ('Name', 'Arithmetic operations with images')
%sum 100 to an image
Iplus100 = imadd(I, 100);
subplot(3, 1, 1);
imshow(Iplus100);
title('Big ben plus 100');

%subtract an image from another one
A = imread('bigben.png');
B = imread('edin_castle.png');
AminusB = imsubtract(A,B);
subplot(3, 1, 2);
imshow(AminusB);
title('Big ben minus Edin castle');

%invert image
inverted = abs(255-A);    
subplot(3, 1, 3);
imshow(inverted);
title('Big ben inverted');

figure('Name', 'tresholding image')
%tresholding image 0.2
for i = 1:size(A, 1)
    for j = 1:size(A, 2)
        for k = 1:size(A, 3)
            if A(i, j, k) > (0.2*255)
                B(i, j, k) = 255;
            else
                B(i, j, k) = 0;
            end
        end
    end
end
subplot(1, 3, 1);
imshow(B);
title('0.2');

%tresholding image 0.5
for i = 1:size(A, 1)
    for j = 1:size(A, 2)
        for k = 1:size(A, 3)
            if A(i, j, k) > (0.5*255)
                B(i, j, k) = 255;
            else
                B(i, j, k) = 0;
            end
        end
    end
end

subplot(1, 3, 2);
imshow(B);
title('0.5');

%tresholding image 0.7
for i = 1:size(A, 1)
    for j = 1:size(A, 2)
        for k = 1:size(A, 3)
            if A(i, j, k) > (0.7*255)
                B(i, j, k) = 255;
            else
                B(i, j, k) = 0;
            end
        end
    end
end

subplot(1, 3, 3);
imshow(B);
title('0.7');

%% Logarithm transform

I = imread('peppers.png'); %load image
Imgray = rgb2gray(I); %convert it to gray scale
Imgraydouble = im2double(Imgray);
figure('Name', 'Logarithm Transform');
% with C = 2
Ic2 = 2.*log10(1+Imgraydouble);

subplot(2, 3, 1);
imshow(Ic2);
title('2');

% with C = 5
Ic5 = 5.*log10(1+Imgraydouble);

subplot(2, 3, 2);
imshow(Ic5);
title('5');

% with C = 7
Ic7 = 7.*log10(1+Imgraydouble);

subplot(2, 3, 3);
imshow(Ic7);
title('7');

% with C = 15
Ic15 = 15.*log10(1+Imgraydouble);

subplot(2, 3, 4);
imshow(Ic15);
title('15');

% with C = 20
Ic20 = 20.*log10(1+Imgraydouble);

subplot(2, 3, 5);
imshow(Ic20);
title('20');

subplot(2, 3, 6);
imshow(Imgray);
title('original');

%% Exponencial trnsform

figure('Name', 'Exponential Transform');
% with alpha = 0.3
Ia3 = 5.*(((1+ 0.3).^Imgraydouble)-1);

subplot(1, 3, 1);
imshow(Ia3);
title('0.3');

% with alpha = 0.4
Ia4 = 5*(((1+ 0.4).^Imgraydouble)-1);

subplot(1, 3, 2);
imshow(Ia4);
title('0.4');

% with alpha = 0.6
Ia6 = 5.*(((1+ 0.6).^Imgraydouble)-1);

subplot(1, 3, 3);
imshow(Ia6);
title('0.6');

%% Power Transformation law
figure('Name', 'Power Law Transform');
% with gamma = 0.5
Ip05 = 2.*(Imgraydouble.^0.5);

subplot(1, 3, 1);
imshow(Ip05);
title('0.5');

% with gamma = 1.5
Ip15 = 2.*(Imgraydouble.^1.5);

subplot(1, 3, 2);
imshow(Ip15);
title('1.5');

% with gamma = 3
Ip30 = 2.*(Imgraydouble.^3);


subplot(1, 3, 3);
imshow(Ip30);
title('3.0');

%% Histogram Stretching
low = 0.05;
high = 0.95; 
scretched = imadjust(Imgray, [low;high], [0;1]);

figure('Name', 'Stretched Image and Histogram')
subplot(1, 2, 1);
imshow(scretched);
title('Stretched Image');

subplot(1, 2, 2);
imhist(scretched);
title('Histogram');

%% Image Equalization

equalized = histeq(Imgray);
figure('Name','Histogram Equalization')
subplot(2, 1, 1);
imshowpair(Imgray, equalized, 'montage');
axis off;
subplot(2, 2, 3);
imhist(Imgray);
title('Original');
subplot(2, 2, 4);
imhist(equalized);
title('Equalized');

%% Histogram matching

matched = histeq(Imgray, (0:255));
figure('Name','Histogram Matching')
subplot(2, 1, 1);
imshowpair(Imgray, matched, 'montage');
axis off;
subplot(2, 2, 3);
imhist(Imgray);
title('Original');
subplot(2, 2, 4);
imhist(matched);
title('Matched');
%% Equalization of a RGB Image
hsvI = rgb2hsv(I);
hsvEqualized = hsvI;
hsvEqualized(:, :, 3) = histeq(hsvI(:,:,3));
rgbEqualized = hsv2rgb(hsvEqualized);
figure('Name','RGB Equalization')
imshowpair(I, rgbEqualized, 'montage');
axis off;

%% Add salt pepper and gausian noise

salt_peper = imnoise(I, 'salt & pepper');



%% Mean filter
hmean = ones(3, 3, 3)/9;
imageMean = imfilter(I, hmean);
