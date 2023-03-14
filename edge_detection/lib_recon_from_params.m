function y = lib_recon_from_params(x, signal_characteristics, param_method)

% Function that reconstructs the signal from the parameters
%
% INPUT
% x - parameters coefficients
% N - length of reconstructing signal
% P - polynomials
%
% OUTPUT:
% y - reconstructed signal from input parameters
% ------------------------------------------------------------------------

N=signal_characteristics.N;
M=signal_characteristics.M;


P=param_method.P;
y=reshape(sum(x.*P,2),M,N); %synthesize the estimated signal





