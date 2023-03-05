%% solving using gradient FBB-PD
function [batch,batch_cnt]=lib_FBB_PD(param_method, param_solver, y , N, degree, batch,batch_cnt)

% Function which minimazes f1(x) + f2(Lx) by FBB-PD (FBB-PD) algorithm.
% This function calls the Unlocbox.
%
% INPUT:
% param_method - structure, which includes parameters of functions f1 and f2.
% param_solver - structure, which includes parameters used in CP
% y - observed signal
% N - length of signal y
% degree - polynom degree
% batch - structure of results
% batch cnt - number of results in batch
%
% OUTPUT:
% batch - structure of results
% batch cnt - number of results in batch

% -------------------------------------------------------------------------
%initialization point for iterations
start_point = zeros(degree*N,1);

batch_cnt = batch_cnt + 1; %number of results in batch increase

% set input parameters of function f1 and f2
A=param_method.A_mat;
vector_tau=param_method.tau;
% set weights of differences - Default: no weights = vector of ones.
if ~isfield(param_method, 'weights') 
    param_method.weights=ones((N-1),degree); 
end
vector_weights=param_method.weights;

% compute L norm
norm_L=4*sum((vector_tau.*max(vector_weights)).^2);

% creation of matrix of tau and weights for next computation 
matrix_of_tau=ones((N-1)*degree,1);
matrix_of_weights=ones((N-1)*degree,1);

for i=1:degree
 matrix_of_tau(1+(N-1)*(i-1):(N-1)*i)=matrix_of_tau(1+(N-1)*(i-1):(N-1)*i).*vector_tau(i);
 matrix_of_weights(1+(N-1)*(i-1):(N-1)*i)=matrix_of_weights(1+(N-1)*(i-1):(N-1)*i).*vector_weights(:,i);
end

% setting data term
fdata.eval = @(x) 0.5*norm(A*x-y)^2;
fdata.grad = @(x) A'*(A*x-y);
fdata.beta = param_method.beta;


% Setting term with linear operator L (function f2)
% Setting L
L = @ (x) matrix_of_tau.*matrix_of_weights.*(reshape((diff(reshape(x,N,degree))),(N-1)*degree,1));
Lt = @(x) reshape(diff([zeros(1,degree);reshape((-matrix_of_tau.*(matrix_of_weights.*x)),N-1,degree);zeros(1,degree)]),degree*N,1);
flinear.L = L; 
flinear.Lt = Lt;
flinear.norm_L = norm_L; 

% Group soft thresholding (without any weight, as they are inside L)
flinear.prox = @(x,T) reshape(prox_l21(reshape(x,N-1,degree),T), degree*(N-1),1);
flinear.eval = @(x) norm(sqrt(sum(reshape(x,N,degree).^2,2)),1);


% setting dummy, null function as the third, non-smooth term
f3.eval= @(x) 0;
f3.prox = @(x,T) x;

%Solving
tic
[sol, info] = fb_based_primal_dual(start_point, fdata, flinear, f3, param_solver);
time_end_FBB_PD=toc;

% Reconstruction of observed parameters from CP algorithm
y_recon = lib_recon_from_params(sol, N, degree, param_method); 

% save into batch
batch{batch_cnt}.algorithm = 'FBD';
batch{batch_cnt}.fdata_smooth = fdata;
batch{batch_cnt}.f_nonsmooth = flinear;
batch{batch_cnt}.sol = sol;
batch{batch_cnt}.info = info;
batch{batch_cnt}.info.objective = info.objective(~isnan(info.objective));  %just reliable values
batch{batch_cnt}.y_recon = y_recon;
batch{batch_cnt}.final_objective = fdata.eval(sol) + flinear.eval(sol);
batch{batch_cnt}.elapsed_time = time_end_FBB_PD;