%% DEMO - re-weighted and non-weighted Condat algorithm
% This script performs the synthetic experiment with greater noise 
%
% What needs to be set:
% It is necessary to use the UnlocBox toolbox
% signal_type - 1 - linear, 2 - polynomial (degree=3)
%
%--------------------------------------------------------------------------
% loaded signal = signal_characteristic - a structrure which includes
%   - N - length of signal
%   - degree - polynomial degree
%   - y - noisy signal
%   - signal_clean - signal without noise
%   - sigma - noise standard deviation
%   - coeff - parameterization coefficients of signal without noise
%   - jumps - number of jumps in the signal
% parameters of solver 
%   - MAX_ITER - maximum number of iterations
% parameters for methods (param_method)
%   - A - representation matrix
%   - weights of differences - Default: vector of ones
%   - tau_x0, tau_x1 and tau_x2
% parameters for breakpoint detection (param_method)
%   - threshold - for breakpoint detection
%   - window_length - length of window in median filter
%   - p - type of norm. Default p=2 -> Euclidean norm.
% %   - proximity - smallest proximity of two breakpoints. Default
%   proximity=2
%--------------------------------------------------------------------------
 
clc
clear variables;
close all;

%% choose type of the signal
% 1 - linear, 2 - polynomial

signal_type=1;

%% Load signal 

if signal_type==1
    load ('test_data\test_signal_4_jmp_6_SNR_25_17.mat')
elseif signal_type==2
    load('test_data\test_signal_13_jmp_7_SNR_25_20.mat')
end

 
% properties of the loaded signal
degree=signal_characteristics.degree;
N=signal_characteristics.N; %length of the signal
y=signal_characteristics.y;
signal_clean=signal_characteristics.signal_clean;
sigma=signal_characteristics.sigma;
njumps=signal_characteristics.jumps;

parameterization_coeff=signal_characteristics.coeff;

%% use of 3-polynomial bases for linear signal

% if degree==2
%     degree=3;
%     signal_characteristics.degree=3;
%     parameterization_coeff=[signal_characteristics.coeff, zeros(N,1)];
% end

%%

% Compute SNR of observed signal y
y_error=y-signal_clean;
SNR_dB_orig = 20*log10(norm(signal_clean)/norm(y_error));

%% Plot signal and its parameterization
handle_signal = figure('Name','Signal','NumberTitle','off');
plot(1:N,signal_clean);
hold on;
plot(1:N,y,'r');
title(['SNR = ', num2str(SNR_dB_orig), ' dB'])
legend('clean', 'observed');

if degree==2
    handle_params = figure('Name','Parameterization coefficients','NumberTitle','off');
    plot(1:N,parameterization_coeff(:,1),'bx')
    hold on
    plot(1:N,parameterization_coeff(:,2),'ro')
    legend('x_0', 'x_1');

elseif degree==3
    
    handle_params = figure('Name','Parameterization coefficients','NumberTitle','off');
    plot(1:N,parameterization_coeff(:,1),'bx')
    hold on
    plot(1:N,parameterization_coeff(:,2),'ro')
    hold on
    plot(1:N,parameterization_coeff(:,3),'g*')
    legend('x_0', 'x_1', 'x_2');
    
end

%% Parameters

% Representation matrix
A=zeros(N,N*degree);
A2=zeros(N,degree);
for i=1:degree
    d=(1:N)/N;
    d=d.^(i-1);
    D=diag(d);
    A(:,1+(i-1)*N:i*N)=D;
    A2(:,i)=d';
end
norm_A=norm(A);
A = sparse(A);  %making it a bit faster with no effort
 
param_method.A=A2; 
param_method.A_mat=A; 
param_method.norm_A = norm_A;
param_method.epsilon_weights=1;

param_solver.tol = 1e-6;    
param_solver.debug_mode = 1;
param_solver.verbose = 0; % display parameter


%% Parameters setting

param_method.epsilon_ball=norm(y-sum(A2.*parameterization_coeff,2));
param_method.epsilon_ball_2=norm(y-signal_clean);

% structure for saving results of used methods
batch = cell(0,0);
batch_cnt = 0;

%% Solve using re-weighted algorithm
% Parameters of solver
MAX_ITER=1000; 
MAX_ITER_Condat=3; 
param_solver.maxit_condat=MAX_ITER_Condat;
param_solver.maxit = MAX_ITER;       % maximum number of iterations

param_method.tau=[5*sigma,3.5*sigma,2.5*sigma]; 

[batch,batch_cnt] = lib_reweighted_Condat_hh(param_method, param_solver, y , N, degree, batch, batch_cnt); 
% recomputed SNR of achieved signal
y_error = batch{batch_cnt}.y_recon-signal_clean;
SNR_dB_Cr_hh = 20*log10(norm(signal_clean)/norm(y_error));


%% Solve using non-weighted algorithm
% Parameters of solver
MAX_ITER=3000;
MAX_ITER_Condat=1;
param_solver.maxit_condat=MAX_ITER_Condat;
param_solver.maxit = MAX_ITER;       % maximum number of iterations

param_method.tau=[5*sigma,3.5*sigma,2.5*sigma];
[batch,batch_cnt] = lib_nonweighted_Condat_hh(param_method, param_solver, y , N, degree, batch, batch_cnt); % recomputed SNR of achieved signal
y_error = batch{batch_cnt}.y_recon-signal_clean;
SNR_dB_C_hh = 20*log10(norm(signal_clean)/norm(y_error));


%% Plot results
colors = 'gcmyrk'; %colors for succesive plots
letters = 'xo*.';

% Display signal
figure(handle_signal);
for cnt = 1:batch_cnt
    plot(batch{cnt}.y_recon, colors(cnt));
end
legend('clean','observed', 'rCA_{hh}', 'nonCA_{hh}');

%% Breakpoint detection and signal denoising
% set parameters
param_method.threshold=0.2;
param_method.window_length=5;  %window length for the median filter
% breakpoints detection and signal denoising using least squares
[batch, batch_cnt]=lib_signal_denoising(y, batch, batch_cnt, N, degree, param_method);
disp(['Breakpoints in the model: ', num2str(find(diff(parameterization_coeff(:,1))~=0)')])

 
figure('Name','Objective functional during iterations','NumberTitle','off'); 
colors = 'brgmcyr';
for cnt = 1:batch_cnt
    plot(batch{cnt}.info.objective, colors(cnt));
        hold on
    
end
axis tight
xlabel('Iterations')
ylabel('Values of objective function')
legend( 'Objective function rCA_{hh}', 'Objective function nonCA_{hh}')

%% statictics evaluation
disp(['Statistics: SNR of observed signal: ', num2str(SNR_dB_orig)])

for cnt = 1:batch_cnt
    batch=lib_evaluation(batch, cnt, signal_characteristics, param_method);
end
