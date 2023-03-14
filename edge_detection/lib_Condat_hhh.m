function [batch,batch_cnt]=lib_Condat_hhh(param_method, param_solver, signal_characteristics, batch, batch_cnt)
N = signal_characteristics.N;
M = signal_characteristics.M;
degree = signal_characteristics.degree;

lambda = param_method.lambda;

start_point = zeros(N*M,degree*degree);
batch_cnt = batch_cnt + 1;

% set input parameters of function f1 and f2 and f3
param.verbose = 3;
param.y = signal_characteristics.signal_noisy;
param.epsilon = param_method.epsilon_ball;

%% f1 - data term 
P = param_method.P;
f1.eval = @(x) eps;
f1.prox = @(x,T) proj_ball(x,param);
f1.L = @(x) sum(P.*x,2);
f1.Lt = @(x) P.*repmat(x,1,degree*degree);
f1.norm_L = param_method.norm_P; %squared norm

%% f2 vertical 
f2.L = @ (x) (reshape((diff(reshape(x,M,N,degree*degree))),(M-1)*N,degree*degree));
f2.Lt = @(x) reshape(diff([zeros(1,N,degree*degree);reshape((-x),M-1,N,degree*degree);zeros(1,N,degree*degree)]),M*N,degree*degree);

f2.norm_L = 4*(degree*degree); %squared norm

f2.prox = @(x,T) prox_l21(x,T, param_solver);
f2.prox_nov = @(x,T) soft_thresh_group(x,T);
f2.eval = @(x) norm(sqrt(sum(x.^2,2)),1);

%% f3 horizontal
f3.L = @ (x) (reshape((diff(reshape(x,M,N,degree*degree),1,2)),M*(N-1),degree*degree));
f3.Lt = @(x) reshape(diff([zeros(M,1,degree*degree),reshape((-x),M,N-1,degree*degree),zeros(M,1,degree*degree)],1,2),M*N,degree*degree);

f3.norm_L=4*degree*degree; %squared norm

f3.prox = @(x,T) prox_l21(x,lambda*T, param_solver);
f3.prox_nov = @(x,T) soft_thresh_group(x,lambda*T);
f3.eval = @(x) norm(sqrt(sum(x.^2,2)),1);

%% Solving
[sol2, info2] = condat_alg_1_hhh(start_point, f1, f2, f3, param_solver);

y_recon = lib_recon_from_params(sol2.primal,signal_characteristics, param_method); 
 
% save into batch
batch{batch_cnt}.algorithm = 'Condat_hhh';
batch{batch_cnt}.sol = sol2.primal;
batch{batch_cnt}.info = info2;
batch{batch_cnt}.signal_recon = y_recon;
batch{batch_cnt}.final_objective = f2.eval(f2.L(sol2.primal)) + f1.eval(f1.L(sol2.primal))+f3.eval(f3.L(sol2.primal)) ;
