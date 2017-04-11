%% Informatica Medica Trabralho 2
 % Processamento de Imagens
 % Pedro Henrique Fornari
 % 13104320

%% Load Image
I = imread('bigben.png');
Igray = rgb2gray(I);
figure ('Name', 'Original Image');
imshow(I);

%% Apply bit plane slicing
B = cell(8);

figure('Name', 'Bit Panel Slice')
for i = 1:8
    B{i}=zeros(size(Igray));
    B{i}=bitset(B{i},i,bitget(Igray,i));
    B{i}=uint8(B{i});
    subplot(2, 4, i);
    imshow(B{i}*255);
    str = sprintf('Panel Bit %d', i);
    title(str);
end

%% Get RGB channels

Rchannel = I(:, :, 1);
Gchannel = I(:, :, 2);
Bchannel = I(:, :, 3);

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

alfa = 0.2989; beta = 0.5870; gamma = 0.1140;
for i = 1:size(I, 1)
    for j = 1:size(I, 2)
        Imgray(i, j) = alfa*I(i, j, 1) + beta*I(i, j, 2) + gamma*I(i, j, 3);
    end
end

figure('Name', 'Manual conversion from rgb to gray scale')
imshow(Imgray);

