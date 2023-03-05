function [batch,batch_cnt]=lib_DR_poly(param_method, param_solver, y , N, degree, batch, batch_cnt)

% Function which minimazes f1(x) + f2(x) by DR algorithm.
% This function calls the Unlocbox.
%
% INPUT:
% param_method - structure, which includes parameters of functions f1 and f2.
% param_solver - structure
% y - observed signal
% N - length of signal y
% degree - polynom degree
% batch - structure of results
% batch cnt - number of results in batch
%
% OUTPUT:
% batch - structure of results
% batch cnt - number of results in batch

batch_cnt = batch_cnt + 1;
start_point = zeros(degree*N,1);

% setting of input parameters
A = param_method.A_mat;
vector_tau=param_method.tau;
param_tv1d = param_solver.param_tv1d;

% data term
fdata.eval = @(x) 0.5*norm(A*x-y)^2;
O=matrix_inversion(A,1);
O2=1*A'*y;
fdata.prox = @(x,T) O*(O2+x);

% setting f1 (nonsmooth)
if degree==2
    f_nonsmooth.prox = @(x,T) [ prox_tv1d(x(1:N),T*vector_tau(1),param_tv1d); prox_tv1d(x(N+1:end),T*vector_tau(2),param_tv1d) ];
elseif degree==3
    f_nonsmooth.prox = @(x,T) [ prox_tv1d(x(1:N),T*vector_tau(1),param_tv1d); prox_tv1d(x(N+1:2*N),T*vector_tau(2),param_tv1d); ...
        prox_tv1d(x(2*N+1:end),T*vector_tau(3),param_tv1d)];
end

f_nonsmooth.eval = @(x) sum(abs(diff(reshape(x,N,degree))))*vector_tau';

%solve using Unlocbox
tic
[sol, info] = douglas_rachford(start_point, f_nonsmooth, fdata, param_solver);
time_end_DR=toc;

y_recon = lib_recon_from_params(sol, N, degree, param_method); %synthesize the estimated signal

% save into batch
batch{batch_cnt}.algorithm = 'DR';
batch{batch_cnt}.fdata = fdata;
batch{batch_cnt}.f_nonsmooth = f_nonsmooth;
batch{batch_cnt}.sol = sol;
batch{batch_cnt}.info = info;
batch{batch_cnt}.info.objective = info.objective(~isnan(info.objective));  %just reliable values
batch{batch_cnt}.y_recon = y_recon;
batch{batch_cnt}.final_objective = fdata.eval(sol) + f_nonsmooth.eval(sol);
batch{batch_cnt}.TVx0 = norm(diff(sol(1:N),1));
batch{batch_cnt}.TVx1 = norm(diff(sol(N+1:end),1));
batch{batch_cnt}.elapsed_time = time_end_DR;
