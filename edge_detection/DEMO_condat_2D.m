% DEMO - edge detection in images

% This script performs the experiment with one loaded image 
% Used algorithm is Condat
%
% What needs to be set:
% It is necessary to use the images from training subset of BSDS300
% database

% need to be set
% - sigma - which defines added noise - 0 for clean images
% - SNR - which defines added noise
% - degree - polynomial degree (K+1)
% - delta - model error
% - type_of_polynomials - type of used polynomials: 1 - standard
% - param_method.lambda - lambda
% - param_solver.maxit - maximum number of iterations 
%
%--------------------------------------------------------------------------
% Other parameters and variables
% - signal_characteristic - a structrure which includes
%   - M - size of image
%   - N - size of image
%   - degree - polynomial degree
%   - signal_noisy - noisy signal
%   - signal_clean - signal without noise
%   - sigma - noise standard deviation
%   - diff - image of differences
%   - signal_recon - reconstructed image
% parameters of solver 
%   - MAX_ITER - maximum number of iterations 
% parameters for edge detection
%   - threshold - for edge detection

%-----------------------------------------------------

clear variables
close all

% ----------------------------------------------------
%Input data:
imgDir = 'data/images';

%Output folders
outDir_gray = 'data/outputs_gray';
outDir_RGB = 'data/outputs_RGB';
outDir_gray_noisy = 'data/outputs_noisy';

%-----------------------------------

degree = 3; % K - polynomial degree
type_of_polynomials = 1; %standard basis
SNR = 20; % SNR - make noisy grayscale image

param_method.lambda = 1; 
param_solver.maxit = 1000;       % maximum number of iterations


iids = dir(fullfile(imgDir,'*.jpg'));

% take first image
name_of_image=iids(1).name;

%% grayscale images
% Experiment 1
delta = 6250; 
sigma = 0;
[signal_characteristics_gray]=lib_add_noise_to_image(sigma, [imgDir, '/', name_of_image]);

% compute edges
param_method.epsilon_ball = delta;
signal_characteristics_gray.degree = degree;
[signal_characteristics_gray] = edge_detect_2D(param_solver, param_method, signal_characteristics_gray, type_of_polynomials);
signal_characteristics_gray.delta=delta;

results_gray=signal_characteristics_gray;

% norm differences
diff = abs(results_gray.diff);
results_gray.diff_norm = (diff-min(min(diff)))./(max(max(diff))- min(min(diff)));
% thinned differences
threshold = 0.125;
results_gray.diff_threshold_thin = double(bwmorph((results_gray.diff_norm>=threshold), 'thin', inf));

% save results
save ([outDir_gray, '/', name_of_image(1:end-4), '.mat'], 'results_gray')

disp('GRAY done')
%% RGB
% Experiment 2
delta = 7500; 
sigma = 0;
[signal_characteristics_RGB]=lib_add_noise_to_image_rgb(sigma, [imgDir, '/', name_of_image]);

% compute edges
param_method.epsilon_ball = delta;
signal_characteristics_RGB.degree = degree;
[signal_characteristics_RGB] = edge_detect_2D_rgb(param_solver, param_method, signal_characteristics_RGB, type_of_polynomials);
signal_characteristics_RGB.delta=delta;

results_RGB=signal_characteristics_RGB;




% norm differences
clear diff
diff = abs(results_RGB.diff);
diff2=max(diff,[],3); % merging variant 1 - maximum value of R,G,B
results_RGB.diff_norm = (diff2-min(min(diff2)))./(max(max(diff2))- min(min(diff2)));
% thinned differences
threshold = 0.125;
results_RGB.diff_threshold_thin = double(bwmorph((results_RGB.diff_norm>=threshold), 'thin', inf));

% save results
save ([outDir_RGB, '/', name_of_image(1:end-4), '.mat'], 'results_RGB')

disp('RGB done')

%% GRAY NOISY
%Experiment 3
delta = 6250; 
[signal_characteristics_noisy]=lib_add_noise_to_image_for_testing(SNR, [imgDir, '/', name_of_image]);

% compute edges
param_method.epsilon_ball = delta;
signal_characteristics_noisy.degree = degree;
[signal_characteristics_noisy] = edge_detect_2D(param_solver, param_method, signal_characteristics_noisy, type_of_polynomials);

results_noisy=signal_characteristics_noisy;

% norm differences
clear diff
diff = abs(results_noisy.diff);
results_noisy.diff_norm = (diff-min(min(diff)))./(max(max(diff))- min(min(diff)));
% thinned differences
threshold = 0.125;
results_noisy.diff_threshold_thin = double(bwmorph((results_noisy.diff_norm>=threshold), 'thin', inf));

% save results
save ([outDir_gray_noisy, '/', name_of_image(1:end-4), '.mat'], 'results_noisy')

disp('noisy done')

%% plot results

% input image
image_clean_rgb = imread([imgDir, '/', name_of_image]);
image_clean_gray = rgb2gray(image_clean_rgb);
image_gray_noisy = (reshape(signal_characteristics_noisy.signal_noisy,signal_characteristics_noisy.M,signal_characteristics_noisy.N));

% image recon
image_recon_rgb = uint8(results_RGB.signal_recon);
image_recon_gray = uint8(results_gray.signal_recon);
image_recon_noisy = uint8(results_noisy.signal_recon);

% differences
diff_norm_gray = results_gray.diff_norm;
diff_norm_rgb = results_RGB.diff_norm;
diff_norm_noisy = results_noisy.diff_norm;

% thinned differences
diff_threshold_thin_gray = results_gray.diff_threshold_thin;
diff_threshold_thin_rgb = results_RGB.diff_threshold_thin;
diff_threshold_thin_noisy = results_noisy.diff_threshold_thin;

figure ('Name', 'Results - RGB + GRAY + LAB')
subplot(3,4,1)
imagesc(image_clean_gray)
colormap(gray)

subplot(3,4,2)
imagesc(image_recon_gray)
colormap(gray)

subplot(3,4,3)
imagesc(diff_norm_gray)
colormap(gray)

subplot(3,4,4)
imagesc(diff_threshold_thin_gray)
colormap(gray)

subplot(3,4,5)
imagesc(image_clean_rgb)
colormap(gray)

subplot(3,4,6)
imagesc(image_recon_rgb)

subplot(3,4,7)
imagesc(diff_norm_rgb)
colormap(gray)

subplot(3,4,8)
imagesc(diff_threshold_thin_rgb)
colormap(gray)

subplot(3,4,9)
imagesc(image_gray_noisy)
colormap(gray)

subplot(3,4,10)
imagesc(image_recon_noisy)
colormap(gray)

subplot(3,4,11)
imagesc(diff_norm_noisy)
colormap(gray)

subplot(3,4,12)
imagesc(diff_threshold_thin_noisy)
colormap(gray)

