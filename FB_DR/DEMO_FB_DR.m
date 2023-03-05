%% DEMO - FB and DR algorithm
%
% This script performs the experiment with one loaded signal 
%
% What needs to be set:
% It is necessary to use the UnlocBox toolbox
% signal_type - 1 - linear, 2 - polynomial (degree=3)
%
%--------------------------------------------------------------------------
% Other parameters and variables
% loaded signal = signal_characteristic - a structrure which includes
%   - N - length of signal
%   - degree - polynomial degree
%   - y - noisy signal
%   - signal_clean - signal without noise
%   - sigma - noise standard deviation
%   - coeff - parameterization coefficients of signal without noise
% parameters of solver 
%   - MAX_ITER - maximum number of iterations
% parameters for methods (param_method)
%   - A - representation matrix
%   - tau 
% parameters for breakpoint detection (param_method)
%   - threshold - for breakpoint detection
%   - window_length - length of window in median filter
%   - p - type of norm. Default p=2 -> Euclidean norm.
%   - proximity - smallest proximity of two breakpoints.
%     Default proximity=2

%% 
clc
clear variables;
close all;

%% choose type of the signal
% 1 - linear, 2 - polynomial

signal_type=2;

%% Load signal 

if signal_type==1
    load ('test_data\test_signal_4_jmp_6_SNR_25_17.mat')
elseif signal_type==2
    load('test_data\test_signal_13_jmp_7_SNR_25_20.mat')
end

% properties of the loaded signal
degree=signal_characteristics.degree;
y=signal_characteristics.y;
signal_clean=signal_characteristics.signal_clean;
parameterization_coeff=signal_characteristics.coeff;
sigma=signal_characteristics.sigma;
N=signal_characteristics.N;

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

%% Plot signal and its parameterization coefficients
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


%% parameters from the minimization formulation and the algorithm
MAX_ITER = 3000;

param_solver.verbose = 0;    % display parameter
param_solver.maxit = MAX_ITER;     % maximum number of iterations
param_solver.tol = 1e-6;        % tolerance to stop iterating (by relative error)
param_solver.debug_mode = 1;
% param_solver.method = 'ISTA';  % desired method for solving the problem. By default, it's 'FISTA' in FB
param_tv1d.use_fast = 'true';  %use fast Condat's algorithm for TV prox
param_solver.param_tv1d.use_fast = 'true';  %use fast Condat's algorithm for TV prox


% Representation matrix
d=(1:N)/N;
A2=zeros(N,degree);

for i = 1:degree
    d2=d.^(i-1);
    D = diag(d2);
    A(:,1+(i-1)*N:i*N) = D;
    A2(:,i)=d2';
end

beta = norm(A)^2; % Lipschitz constant %norm(A)^2;
A = sparse(A);  %making it a bit faster with no effort

param_method.beta=beta;
param_method.A=A2;
param_method.A_mat=A;


%% Solving

%prepare a cell array of methods (and results)
batch = cell(0,0);
batch_cnt = 0;

%% Solving using gradient (FB)

if degree==2
    param_method.tau = [sigma*5, sigma*3.5];
elseif degree==3
    param_method.tau = [sigma*5, sigma*3.5, sigma*2.5];
end

[batch,batch_cnt]=lib_FB_poly(param_method, param_solver, y , N, degree, batch, batch_cnt);


% recomputed SNR of achieved signal
y_error=batch{batch_cnt}.y_recon-signal_clean;
SNR_dB_FB = 10*log10(sum(signal_clean.^2)/sum(y_error.^2));

%% Solving not using gradient (DR)
param_solver.gamma = 0.01;   % stepsize

if degree==2
    param_method.tau = [sigma*200, sigma*150];
elseif degree==3
    param_method.tau = [sigma*200, sigma*150, sigma*130];
end

[batch,batch_cnt]=lib_DR_poly(param_method, param_solver, y , N, degree, batch, batch_cnt);

y_error=batch{batch_cnt}.y_recon-signal_clean;
SNR_dB_DR = 10*log10(sum(signal_clean.^2)/sum(y_error.^2));

%% display the signals and parameters

colors = 'gcmy'; %colors for succesive plots

%signal
figure(handle_signal);
for cnt = 1:batch_cnt
    plot(batch{cnt}.y_recon, colors(cnt));
end
legend('clean','observed','FB', 'DR');

figure ('Name', 'Objective function', 'NumberTitle', 'off')
for cnt = 1:batch_cnt
    plot(batch{cnt}.info.objective, colors(cnt));
    hold on
end
set(gca,'YScale','log')
axis tight
title('Objective functional during iterations')
legend('objective function FB', 'objective function DR')


%% Breakpoint detection and signal denoising
% set parameters
param_method.threshold = 0.2;
param_method.window_length = 5;  %window length for the median filter
param_method.p = 2;  %window length for the median filter

% breakpoints detection and signal denoising using least squares
[batch, batch_cnt] = lib_signal_denoising(y, batch, batch_cnt, N, degree, param_method);
disp(['Breakpoints in the model: ', num2str(find(diff(parameterization_coeff(:,1))~=0)')])

%% ploting of signal reconstruction
y_error=batch{1}.signal_recon-signal_clean;
SNR_dB_FB_recon = 10*log10(sum(signal_clean.^2)/sum(y_error.^2));

y_error=batch{2}.signal_recon-signal_clean;
SNR_dB_DR_recon = 10*log10(sum(signal_clean.^2)/sum(y_error.^2));

figure ('Name', 'Comparison - reconstructed signal', 'NumberTitle', 'off')
title('Comparison - reconstructed signal')
plot(y,'r')
hold on
plot(signal_clean)
hold on
for cnt = 1:batch_cnt
    plot(batch{cnt}.signal_recon, colors(cnt));
    hold on
end
legend('signal noisy', 'signal clean', 'signal recon FB', 'signal recon DR');

%% statictics evaluation
disp(['Statistics: SNR of observed signal: ', num2str(SNR_dB_orig)])
    

for cnt = 1:batch_cnt
    batch=lib_evaluation(batch, cnt, signal_characteristics, param_method);
end

