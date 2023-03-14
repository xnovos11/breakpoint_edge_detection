function [signal_characteristics]=edge_detect_2D_rgb(param_solver, param_method, signal_characteristics, type_of_polynomials)
%% This script performs the experiment with noisy image 
% loaded signal = signal_characteristics - a structrure which includes
%   - N - size of image
%   - M - size of image
%   - degree - polynomial degree
%   - signal_clean - signal without noise
%   - signal_noise - signal with noise
%   - sigma - noise standard deviation
% parameters of solver (param_solver)
%   - maxit - maximum number of iterations
% parameters for methods (param_method)
%   - p - representation matrix with polynomials
%   - tau - vector of taus. Default vector of Ones
%   - norm_P - norm of representation matrix P
% parameters for edges detection (param_method)
%   - window_length - length of window in median filter. 
% parameters for polynomials generation (param_generation_polynomials)
%   - plot_polynomials: 0 - do not plot, 1 - plot, DEFAULT: 0
%   - type: 1 - standard, 4 - orthogonal
%--------------------------------------------------------------------------

param_generation_polynomials.plot_polynomials = 1;  
param_generation_polynomials.type = type_of_polynomials;
   
%% Setting of parameters of method
param_method.window_length = 5;  %window length for the median filter

%% Parameters setting

% Representation matrix
[P] = lib_generate_polynomials_2D_ortho(signal_characteristics,param_generation_polynomials);

param_method.P = P; % Matrix of polynomials
param_method.norm_P = max(sum(P.^2,2)); % square norm of P
param_method.p = 2; 

%% structure for saving results of used methods
batch = cell(0,0);
% batch_cnt = 0;
 degree = signal_characteristics.degree;
 N = signal_characteristics.N;
 M = signal_characteristics.M;

%% Solve using CP algorithm 
figure ('Name', 'Obj function')
for batch_cnt=0:2
[batch,batch_cnt] = lib_Condat_hhh_rgb(param_method, param_solver, signal_characteristics, batch, batch_cnt); 
signal_characteristics.signal_recon(:,:,batch_cnt) = batch{1,batch_cnt}.signal_recon; % save reconstructed image


subplot(1,3,batch_cnt)
plot(batch{1,batch_cnt}.info.objective)


%% Edge detection
% Gradients in differneces of parameters

 parameterization_coeff_sol = reshape(batch{1,batch_cnt}.sol,M,N,degree*degree); % reshape given parameters
 
% better contrast for images
 minimum = min(min(min(abs(parameterization_coeff_sol))));
 maximum = max(max(max(abs(parameterization_coeff_sol))));
 parameterization_coeff_sol_auto_contrast = 0 + (abs(parameterization_coeff_sol) - minimum).*((255-0)/(maximum-minimum));
     
 for i=1:degree*degree
     figure ('Name', ['param_coeff_' num2str(i)]); image(parameterization_coeff_sol(:,:,i),'CDataMapping','scaled');
     colormap(gray)
% 
%      imwrite(uint8(parameterization_coeff_sol_auto_contrast(:,:,i)),['fig_param_coeff_' num2str(i) '.png'])
   
     matrix_of_gradients(:,:,i)=abs(imgradient(parameterization_coeff_sol(:,:,i)));
 end
 
 diff_parameterization_coeff = (sum(matrix_of_gradients.^param_method.p,3)).^(1/param_method.p);
 batch{1,batch_cnt}.diff_euclidean = diff_parameterization_coeff;
 signal_characteristics.diff(:,:,batch_cnt) = batch{1,batch_cnt}.diff_euclidean; % save gradients 

end
