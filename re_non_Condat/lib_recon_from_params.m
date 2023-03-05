function y = lib_recon_from_params(x, N, degree, param_method)

% Function that reconstructs the signal from the parameters
%
% INPUT
% x - parameters coefficients
% N - length of reconstructing signal
% degree - polynom degree
%
% OUTPUT:
% y - reconstructed signal from input parameters
% ------------------------------------------------------------------------

A = param_method.A;
y = sum(reshape(x,N,degree).*A,2); %synthesize the estimated signal


