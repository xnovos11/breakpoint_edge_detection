function [signal_characteristics]=lib_add_noise_to_image(sigma, name_of_image)
% This function adds noise to the load image

y = imread(name_of_image);
y=rgb2gray(y);
M = size(y,1); % rows
N = size(y,2); % columns
signal_clean = reshape(double(y), N*M, 1);

% characteristic of image
signal_characteristics.N = N;
signal_characteristics.M = M;
signal_characteristics.signal_clean=signal_clean;


%% add noise
signal_noisy = signal_clean + sigma * randn(N*M,1);  %add noise

signal_characteristics.signal_noisy=signal_noisy;

signal_characteristics.sigma=sigma;

              





