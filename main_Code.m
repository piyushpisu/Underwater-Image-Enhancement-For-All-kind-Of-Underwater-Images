
%% Load the image and split channels. 

rgbImage=double(imread("C:\Users\piyus\OneDrive\Desktop\Project Phase 2 Update\Input Images\sample_5.jpg"))/255;

% Extracting  the individual red, green, and blue color channels because
% we can't white balance a single color channel ( Because a single channel is monochrome, not color).
redChannel = rgbImage(:,:,1);
greenChannel = rgbImage(:,:,2);
blueChannel = rgbImage(:,:,3);

Ir_mean = mean(redChannel, 'all');
Ig_mean = mean(greenChannel, 'all');
Ib_mean = mean(blueChannel, 'all');

%% Color compensation and Correction
% due to a very small mean
% value for the red channel, leading to an overcompensation of
% this channel in locations where red is present (because Gray
% world devides each channel by its mean value). To circumvent
% this issue, we therefore primarily aim to compensate
% for the loss of the red channel.
alpha = .1;
Irc = redChannel + alpha*(Ig_mean - Ir_mean);
alpha = 0; % 0 does not compensates blue channel. 

Ibc = blueChannel + alpha*(Ig_mean - Ib_mean);

%% White Balance

I = cat(3, Irc, greenChannel, Ibc);
I_lin = rgb2lin(I);
percentiles = 25;
illuminant = illumgray(I_lin,percentiles);
I_lin = chromadapt(I_lin,illuminant,'ColorSpace','linear-rgb');
Iwb = lin2rgb(I_lin);

%% Gamma Correction
Igamma = imadjust(Iwb,[],[],2);

%% image sharpening
sigma = 20;
Igauss = Iwb;
N = 25;
for iter=1: N
   Igauss =  imgaussfilt(Igauss,sigma);
   Igauss = min(Iwb, Igauss);
end

gain = .1; %sometimes gain <1 is better. 
Norm = (Iwb-gain*Igauss);

%Normalization
for n = 1:3
   Norm(:,:,n) = histeq(Norm(:,:,n)); 
end
Isharp = (Iwb + Norm)/1.7;

%% weights calculation

% Lapacian contrast weight 
Isharp_lab = rgb2lab(Isharp);
Igamma_lab = rgb2lab(Igamma);

% input1
R1 = double(Isharp_lab(:, :, 1)) / 255;

% calculate laplacian contrast weight
WC1 = sqrt((((Isharp(:,:,1)) - (R1)).^2 +((Isharp(:,:,2)) - (R1)).^2 + ((Isharp(:,:,3)) - (R1)).^2)/3);

% calculate the saliency weight
WS1 = saliency_detection(Isharp);
WS1 = WS1/max(WS1,[],'all');

% calculate the saturation weight
WSAT1 = sqrt(1/3*((Isharp(:,:,1)-R1).^2 + (Isharp(:,:,2)-R1).^2 + (Isharp(:,:,3)-R1).^2));

% input2
R2 = double(Igamma_lab(:, :, 1)) / 255;

% calculate laplacian contrast weight
WC2 = sqrt((((Igamma(:,:,1)) - (R2)).^2 + ((Igamma(:,:,2)) - (R2)).^2 + ((Igamma(:,:,3)) - (R2)).^2)/3);

% calculate the saliency weight
WS2 = saliency_detection(Igamma);
WS2 = WS2/max(WS2,[],'all');

% calculate the saturation weight
WSAT2 = sqrt(1/3*((Igamma(:,:,1)-R1).^2+(Igamma(:,:,2)-R1).^2+(Igamma(:,:,3)-R1).^2));

% calculate the normalized weight
W1 = (WC1 + WS1 + WSAT1+0.1) ./ (WC1 + WS1 + WSAT1 + WC2 + WS2 + WSAT2+0.2);
W2 = (WC2 + WS2 + WSAT2+0.1) ./ (WC1 + WS1 + WSAT1 + WC2 + WS2 + WSAT2+0.2);
 
 
%% Multi scale trasform fusion.
 img1 = Isharp;
 img2 = Igamma;

% calculate the gaussian pyramid
level = 5;
Weight1 = gaussian_pyramid(W1, level);
Weight2 = gaussian_pyramid(W2, level);

% calculate the laplacian pyramid

% input1
R1 = laplacian_pyramid(Isharp(:, :, 1), level);
G1 = laplacian_pyramid(Isharp(:, :, 2), level);
B1 = laplacian_pyramid(Isharp(:, :, 3), level);

% input2
R2 = laplacian_pyramid(Igamma(:, :, 1), level);
G2 = laplacian_pyramid(Igamma(:, :, 2), level);
B2 = laplacian_pyramid(Igamma(:, :, 3), level);

%% Fusion
for k = 1 : level
   Rr{k} = Weight1{k} .* R1{k} + Weight2{k} .* R2{k};
   Rg{k} = Weight1{k} .* G1{k} + Weight2{k} .* G2{k};
   Rb{k} = Weight1{k} .* B1{k} + Weight2{k} .* B2{k};
end

for i = 1 : level
    r_py{i} = Weight1{i} .* R1{i} + Weight2{i} .* R2{i};
    g_py{i} = Weight1{i} .* G1{i} + Weight2{i} .* G2{i};
	b_py{i} = Weight1{i} .* B1{i} + Weight2{i} .* B2{i};
end

for i = level : -1 : 2
    [m, n] = size(r_py{i - 1});
    r_py{i - 1} = r_py{i - 1} + imresize(r_py{i}, [m, n]);
end
R = r_py{1};

for i = level : -1 : 2
    [m, n] = size(g_py{i - 1});
    g_py{i - 1} = g_py{i - 1} + imresize(g_py{i}, [m, n]);
end
G = g_py{1};

for i = level : -1 : 2
    [m, n] = size(b_py{i - 1});
    b_py{i - 1} = b_py{i - 1} + imresize(b_py{i}, [m, n]);
end
B = b_py{1};

% % reconstruct & output
% R = pyramid_reconstruct(Rr);
% G = pyramid_reconstruct(Rg);
% B = pyramid_reconstruct(Rb);

fusion = cat(3, R, G, B);

%% Results

figure('Name','Step I-IV');

subplot(221);
imshow(I);
title('Original Image');

subplot(222);
imshow(Iwb);
title('I. White Balance');

subplot(223);
imshow(Igamma);
title('II. Gamma Corrected');

subplot(224);
imshow(Isharp);
title('III. Sharpened');

figure('Name','Final Comparison');

subplot(121);
imshow(I);
title('Original');

subplot(122);
imshow(fusion);
title('Final Image');