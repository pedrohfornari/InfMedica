% Universidade Federal de Santa Catarina
% Centro Tecnológico
% Departamento de Engenharia Elétrica e Eletrônica
% Disciplina: Introdução à Informática Médica (EEL7307)
% Professora: Christine Fredel Boos
% Alunos: Lucas Pereira Luiz - 13101258
%         Pedro Henrique Kappler Fornari - 13104320
% Semestre letivo: 2017/1
%
% Atividade_02 : Processamento Digital de Imagem
%

clc;
clear;
close all;


%% ----- 1 - BitSlicing -----
info_pure = imfinfo('bigben.png');
img_pure = imread('bigben.png');

img_gray = rgb2gray(img_pure);


figure('Name','Bitslicing','NumberTitle','off')

subplot(2,4,1);
imshow(bitget(img_gray,8)*255);
title('bit 0');

subplot(2,4,2);
imshow(bitget(img_gray,7)*255);
title('bit 1');

subplot(2,4,3);
imshow(bitget(img_gray,6)*255);
title('bit 2');

subplot(2,4,4);
imshow(bitget(img_gray,5)*255);
title('bit 3');

subplot(2,4,5);
imshow(bitget(img_gray,4)*255);
title('bit 4');

subplot(2,4,6);
imshow(bitget(img_gray,3)*255);
title('bit 5');

subplot(2,4,7);
imshow(bitget(img_gray,2)*255);
title('bit 6');

subplot(2,4,8);
imshow(bitget(img_gray,1)*255);
title('bit 7');


%% ----- 2 - RGB Splitting -----

figure('Name','RGB Splittig','NumberTitle','off')

subplot(2,2,1);
imshow(img_pure);
title('Original');

subplot(2,2,2);
imshow(img_pure(:,:,1));
title('Red Channel');

subplot(2,2,3);
imshow(img_pure(:,:,2));
title('Green Channel');

subplot(2,2,4);
imshow(img_pure(:,:,3));
title('Blue Channel');


%% ----- 3 - Manual Grayscale Conversion -----

img_gray_manual = 0.2989*img_pure(:,:,1) + 0.5870*img_pure(:,:,2) + 0.1140*img_pure(:,:,3);


figure('Name','Manual Grayscale Conversion','NumberTitle','off')

subplot(2,1,1);
imshow(img_pure);
title('Original');

subplot(2,2,3);
imshow(img_gray);
title('Auto Grayscale');

subplot(2,2,4);
imshow(img_gray_manual);
title('Manual Grayscale');


%% ----- 4 - HSV Splitting -----

img_hsv = rgb2hsv(img_pure);


figure('Name','HSV Splittig','NumberTitle','off')

subplot(2,2,1);
imshow(img_pure);
title('Original');

subplot(2,2,2);
img_hsv1 = ones(480,640,3);
img_hsv1(:,:,1) = img_hsv(:,:,1);
imshow(hsv2rgb(img_hsv1));
title('Hue');

subplot(2,2,3);
imshow(img_hsv(:,:,2));
title('Saturation');

subplot(2,2,4);
imshow(img_hsv(:,:,3));
title('Value');


%% ----- 5 - Other Opperations -----

% Add a constant:
figure('Name','Add a constant','NumberTitle','off');
subplot(1,2,1);
imshow(img_pure);
title('Original image');

subplot(1,2,2);
imshow(img_pure + (-77));
title('Adding a constant (-77)');


% Subtract two images:
info_castle = imfinfo('edin_castle.png');
img_castle = imread('edin_castle.png');

figure('Name','Subtract two images','NumberTitle','off');
subplot(2,2,1);
imshow(img_castle);
title('Castle');

subplot(2,2,2);
imshow(img_pure);
title('BigBen');

subplot(2,1,2);
imshow(img_pure - img_castle);
title('BigBen - Castle');


% Image negative
figure('Name','Image negative','NumberTitle','off');
subplot(1,2,1);
imshow(img_pure);
title('Original');

subplot(1,2,2);
imshow(abs(255 - img_pure));
title('Negative');


% Thresholding

% 0.2
img_02 = img_pure;
for i = 1:(size(img_pure,1))
    for j = 1:(size(img_pure,2))
        for k = 1:(size(img_pure,3))
            if img_pure(i,j,k) > round(0.2*255);
                img_02 (i,j,k) = 255;
            else
                img_02 (i,j,k) = 0;
            end
        end
    end
end

% 0.5
img_05 = img_pure;
for i = 1:(size(img_pure,1))
    for j = 1:(size(img_pure,2))
        for k = 1:(size(img_pure,3))
            if img_pure(i,j,k) > round(0.5*255);
                img_05 (i,j,k) = 255;
            else
                img_05 (i,j,k) = 0;
            end
        end
    end
end

% 0.7
img_07 = img_pure;
for i = 1:(size(img_pure,1))
    for j = 1:(size(img_pure,2))
        for k = 1:(size(img_pure,3))
            if img_pure(i,j,k) > round(0.7*255);
                img_07 (i,j,k) = 255;
            else
                img_07 (i,j,k) = 0;
            end
        end
    end
end


figure('Name','Thresholding','NumberTitle','off');
subplot(2,2,1);
imshow(img_pure);
title('Original');

subplot(2,2,2);
imshow(img_02);
title('Threshold = 0.2');

subplot(2,2,3);
imshow(img_05);
title('Threshold = 0.5');

subplot(2,2,4);
imshow(img_07);
title('Threshold = 0.7');

%% ----- 6 - Transforms -----

% Log transform
img_log2 = 2*log10(1+im2double(img_gray));
img_log5 = 5*log10(1+im2double(img_gray));
img_log7 = 7*log10(1+im2double(img_gray));
img_log15 = 15*log10(1+im2double(img_gray));
img_log20 = 20*log10(1+im2double(img_gray));

figure('Name','Log transform','NumberTitle','off');
subplot(2,3,1);
imshow(img_gray);
title('Original');

subplot(2,3,2);
imshow(img_log2);
title('c = 2');

subplot(2,3,3);
imshow(img_log5);
title('c = 5');

subplot(2,3,4);
imshow(img_log7);
title('c = 7');

subplot(2,3,5);
imshow(img_log15);
title('c = 15');

subplot(2,3,6);
imshow(img_log20);
title('c = 20');


% Exponential transform

img_exp03 = 5*(((1.3).^img_gray) -1);
img_exp04 = 5*(((1.4).^img_gray) -1);
img_exp06 = 5*(((1.6).^img_gray )-1);

figure('Name','Exponential transform','NumberTitle','off');
subplot(2,2,1);
imshow(img_gray);
title('Original');

subplot(2,2,2);
imshow(img_exp03);
title('\alpha = 0.3');

subplot(2,2,3);
imshow(img_exp04);
title('\alpha = 0.4');

subplot(2,2,4);
imshow(img_exp06);
title('\alpha = 0.6');


% Power-law transform

img_pl05 = 2*(im2double(img_gray).^0.5);
img_pl15 = 2*(im2double(img_gray).^1.5);
img_pl30 = 2*(im2double(img_gray).^3.0);

figure('Name','Power-law transform','NumberTitle','off');
subplot(2,2,1);
imshow(img_gray);
title('Original');

subplot(2,2,2);
imshow(img_pl05);
title('\gamma = 0.5');

subplot(2,2,3);
imshow(img_pl15);
title('\gamma = 1.5');

subplot(2,2,4);
imshow(img_pl30);
title('\gamma = 3.0');


%% ----- 7 - Histogram Stretching ------

img_gray_stretch = imadjust(img_gray,[0.05; 0.95],[0 1]);

figure('Name','Histogram stretching','NumberTitle','off');
subplot(2,2,1);
imshow(img_gray);
title('Original image (grayscale)');

subplot(2,2,2);
imhist(img_gray);
title('Original image histogram');

subplot(2,2,3);
imshow(img_gray_stretch);
title('Stretched image');

subplot(2,2,4);
imhist(img_gray_stretch);
title('Stretched image histogram');


%% ----- 8 - Image equalization -----

img_gray_eq = histeq(img_gray);

figure('Name','Image equaliization','NumberTitle','off');
subplot(2,2,1);
imshow(img_gray);
title('Original image (grayscale)');

subplot(2,2,2);
imhist(img_gray);
title('Original image histogram');

subplot(2,2,3);
imshow(img_gray_eq);
title('Equalized image');

subplot(2,2,4);
imhist(img_gray_eq);
title('Equalized image histogram');


%% ----- 9 - Histogram Matching -----


img_gray_match = histeq(img_gray, (0:255));

figure('Name','Histogram Matching','NumberTitle','off');
subplot(2,2,1);
imshow(img_gray);
title('Original image (grayscale)');

subplot(2,2,2);
imhist(img_gray);
title('Original image histogram');

subplot(2,2,3);
imshow(img_gray_match);
title('Image with matching histogram');

subplot(2,2,4);
imhist(img_gray_match);
title('Matched image histogram');


%% ----- 10 - RGB equalization -----

img_hsv_eq = img_hsv;
img_hsv_eq(:,:,3) = histeq(img_hsv(:,:,3));
img_rgb_eq = hsv2rgb(img_hsv_eq);

figure('Name','RGB equalization','NumberTitle','off');
subplot(1,2,1);
imshow(img_pure);
title('Original image');

subplot(1,2,2);
imshow(img_rgb_eq);
title('Equalized image');


%% ----- 11 - Filtering -----

img_snp = imnoise(img_gray,'salt & pepper',0.03);
img_gauss = imnoise(img_gray,'gaussian',0,0.02);

av_filter = ones(3,3)/9;
rank_domain = ones(5,5);
gauss_filter = fspecial('gaussian',[5 5],2);


% Average filtering
img_gray_av = imfilter(img_gray,av_filter);
img_snp_av = imfilter(img_snp,av_filter);
img_gauss_av = imfilter(img_gauss,av_filter);

figure('Name','Average filtering','NumberTitle','off');
subplot(2,3,1);
imshow(img_gray);
title('Original image');

subplot(2,3,2);
imshow(img_snp);
title('Salt & pepper noise image');

subplot(2,3,3);
imshow(img_gauss);
title('Gaussian noise image');

subplot(2,3,4);
imshow(img_gray_av);
title('Original image filtered (average)');

subplot(2,3,5);
imshow(img_snp_av);
title('Salt & pepper noise image filtered (average)');

subplot(2,3,6);
imshow(img_gauss_av);
title('Gaussian noise image filtered (average)');

% Median filtering
img_gray_md = medfilt2(img_gray,[3 3]);
img_snp_md = medfilt2(img_snp,[3 3]);
img_gauss_md = medfilt2(img_gauss,[3 3]);

figure('Name','Median filtering','NumberTitle','off');
subplot(2,3,1);
imshow(img_gray);
title('Original image');

subplot(2,3,2);
imshow(img_snp);
title('Salt & pepper noise image');

subplot(2,3,3);
imshow(img_gauss);
title('Gaussian noise image');

subplot(2,3,4);
imshow(img_gray_md);
title('Original image filtered (median)');

subplot(2,3,5);
imshow(img_snp_md);
title('Salt & pepper noise image filtered (median)');

subplot(2,3,6);
imshow(img_gauss_md);
title('Gaussian noise image filtered (median)');

% Rank filtering
img_gray_rf = ordfilt2(img_gray,25,rank_domain);
img_snp_rf = ordfilt2(img_snp,25,rank_domain);
img_gauss_rf = ordfilt2(img_gauss,25,rank_domain);

figure('Name','Rank filtering','NumberTitle','off');
subplot(2,3,1);
imshow(img_gray);
title('Original image');

subplot(2,3,2);
imshow(img_snp);
title('Salt & pepper noise image');

subplot(2,3,3);
imshow(img_gauss);
title('Gaussian noise image');

subplot(2,3,4);
imshow(img_gray_rf);
title('Original image filtered (Rank)');

subplot(2,3,5);
imshow(img_snp_rf);
title('Salt & pepper noise image filtered (Rank)');

subplot(2,3,6);
imshow(img_gauss_rf);
title('Gaussian noise image filtered (Rank)');

% Gaussian filtering
img_gray_gf = imfilter(img_gray,gauss_filter);
img_snp_gf = imfilter(img_snp,gauss_filter);
img_gauss_gf = imfilter(img_gauss,gauss_filter);

figure('Name','Gaussian filtering','NumberTitle','off');
subplot(2,3,1);
imshow(img_gray);
title('Original image');

subplot(2,3,2);
imshow(img_snp);
title('Salt & pepper noise image');

subplot(2,3,3);
imshow(img_gauss);
title('Gaussian noise image');

subplot(2,3,4);
imshow(img_gray_gf);
title('Original image filtered (Gaussian)');

subplot(2,3,5);
imshow(img_snp_gf);
title('Salt & pepper noise image filtered (Gaussian)');

subplot(2,3,6);
imshow(img_gauss_gf);
title('Gaussian noise image filtered (Gaussian)');

%% ----- 12 - Border detection -----

lp_filter = fspecial('laplacian');

img_gray_rb = edge(img_gray,'roberts');
img_gray_pw = edge(img_gray,'prewitt');
img_gray_sb = edge(img_gray,'sobel');
img_gray_lp = imfilter(img_gray,lp_filter,'symmetric');

img_snp_rb = edge(img_snp,'roberts');
img_snp_pw = edge(img_snp,'prewitt');
img_snp_sb = edge(img_snp,'sobel');
img_snp_lp = imfilter(img_snp,lp_filter,'symmetric');

img_gauss_rb = edge(img_gauss,'roberts');
img_gauss_pw = edge(img_gauss,'prewitt');
img_gauss_sb = edge(img_gauss,'sobel');
img_gauss_lp = imfilter(img_gauss,lp_filter,'symmetric');

figure('Name','Border detection filtering','NumberTitle','off');
subplot(3,5,1);
imshow(img_gray);
title('Original');

subplot(3,5,2);
imshow(img_gray_rb);
title('Original (Roberts)');

subplot(3,5,3);
imshow(img_gray_pw);
title('Original (Prewitt)');

subplot(3,5,4);
imshow(img_gray_sb);
title('Original (Sobel)');

subplot(3,5,5);
imshow(img_gray_lp);
title('Original (Laplacian)');

subplot(3,5,6);
imshow(img_snp);
title('Salt & pepper');

subplot(3,5,7);
imshow(img_snp_rb);
title('Salt & pepper (Roberts)');

subplot(3,5,8);
imshow(img_snp_pw);
title('Salt & pepper (Prewitt)');

subplot(3,5,9);
imshow(img_snp_sb);
title('Salt & pepper (Sobel)');

subplot(3,5,10);
imshow(img_snp_lp);
title('Salt & pepper (Laplacian)');

subplot(3,5,11);
imshow(img_gauss);
title('Gaussian');

subplot(3,5,12);
imshow(img_gauss_rb);
title('Gaussian (Roberts)');

subplot(3,5,13);
imshow(img_gauss_pw);
title('Gaussian (Prewitt)');

subplot(3,5,14);
imshow(img_gauss_sb);
title('Gaussian (Sobel)');

subplot(3,5,15);
imshow(img_gauss_lp);
title('Gaussian (Laplacian)');


%% ----- 13 - Border enhancing -----

img_gray_lp_sharp = imsubtract(img_gray,img_gray_lp);

unsharp_filter = fspecial('unsharp');
img_gray_unsharp = imfilter(img_gray,unsharp_filter);

figure('Name','Border enhancing','NumberTitle','off');
subplot(1,3,1);
imshow(img_gray);
title('Original image');

subplot(2,3,2);
imshow(img_gray_lp);
title('Laplacian border');

subplot(2,3,3);
imshow(img_gray_lp_sharp);
title('Laplacian sharpening');

subplot(2,3,5);
imshow(img_gray - img_gray_unsharp);
title('Unsharp mask border');

subplot(2,3,6);
imshow(img_gray_unsharp);
title('Unsharp mask filter');

% 6: exponential certo?