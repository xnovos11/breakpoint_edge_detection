function [signal_characteristics]=lib_add_noise_to_image_rgb(sigma, name_of_image)
% This function adds noise to the load image

y = imread(name_of_image);
M = size(y,1); % rows
N = size(y,2); % columns
signal_clean(:,:,1) = reshape(double(y(:,:,1)), N*M, 1);
signal_clean(:,:,2) = reshape(double(y(:,:,2)), N*M, 1);
signal_clean(:,:,3) = reshape(double(y(:,:,3)), N*M, 1);

% characteristic of image
signal_characteristics.N = N;
signal_characteristics.M = M;
signal_characteristics.signal_clean=signal_clean;

signal_characteristics.sigma=sigma;





