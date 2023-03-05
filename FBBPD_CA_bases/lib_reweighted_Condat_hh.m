function [batch,batch_cnt]=lib_reweighted_Condat_hh(param_method, param_solver, y , N, degree, batch, batch_cnt)
%global N param_solver batch batch_cnt y A signal_clean

start_point = ones(degree*N,1);
batch_cnt = batch_cnt + 1;

% set input parameters of function f1 and f2
A=param_method.A;
vector_tau=param_method.tau;
epsilon=param_method.epsilon_weights;
% set weights of differences - Default: no weights = vector of ones.
if ~isfield(param_method, 'weights') 
    param_method.weights=ones((N-1),degree); 
end
vector_weights=param_method.weights;

% creation of matrix of tau and weights for next computation 
matrix_of_tau=ones((N-1)*degree,1);
matrix_of_weights=ones((N-1)*degree,1);

for i=1:degree
 matrix_of_tau(1+(N-1)*(i-1):(N-1)*i)=matrix_of_tau(1+(N-1)*(i-1):(N-1)*i).*vector_tau(i);
 matrix_of_weights(1+(N-1)*(i-1):(N-1)*i)=matrix_of_weights(1+(N-1)*(i-1):(N-1)*i).*vector_weights(:,i);
end


param.verbose = 3;
param.y=y;
param.epsilon = param_method.epsilon_ball*1.05;


f1.eval = @(x) eps;
f1.prox = @(x,T) proj_ball(x,param);
f1.L = @(x) sum(A.*reshape(x,N,degree),2);
f1.Lt = @(x) reshape(A,N*degree,1).*repmat(x,degree,1);
f1.norm_L=param_method.norm_A; %

% setting term with linear operator L
L = @ (x) matrix_of_tau.*matrix_of_weights.*(reshape((diff(reshape(x,N,degree))),(N-1)*degree,1));
Lt = @(x) reshape(diff([zeros(1,degree);reshape((-matrix_of_tau.*(matrix_of_weights.*x)),N-1,degree);zeros(1,degree)]),degree*N,1);

f2.L = L; 
f2.Lt = Lt;
f2.norm_L=2*max(max(reshape(matrix_of_weights,N-1,degree)))*norm(vector_tau(1:degree));
%group soft thresholding (without any weight, as they are inside L)
f2.prox = @(x,T) reshape(prox_l21(reshape(x,N-1,degree),T), degree*(N-1),1);
f2.eval = @(x) norm(sqrt(sum(reshape(x,N-1,degree).^2,2)),1);

%Solving
obj_f=nan(param_solver.maxit*param_solver.maxit_condat+1,1);
obj_f(1)=f1.eval(f1.L(start_point))+f2.eval(f2.L(start_point));
tic
for ii=1:param_solver.maxit_condat
    [sol2, info2] = condat_alg_1_hh(start_point, f1, f2, param_solver);
    y_recon = lib_recon_from_params(sol2.primal,N, degree, param_method); 
p=2;
    new_weights=(1./(sum(diff(reshape(sol2.primal,N,degree)).^p,2).^(1/p)+epsilon));
    obj_f((ii-1)*param_solver.maxit+1+1:ii*param_solver.maxit+1)=info2.objective;
   
    matrix_of_weights=[];
    for i=1:degree
        matrix_of_weights=[matrix_of_weights;new_weights];
    end

    L = @ (x) matrix_of_tau.*matrix_of_weights.*(reshape((diff(reshape(x,N,degree))),(N-1)*degree,1));
    Lt = @(x) reshape(diff([zeros(1,degree);reshape((-matrix_of_tau.*(matrix_of_weights.*x)),N-1,degree);zeros(1,degree)]),degree*N,1);
    f2.L = L; 
    f2.Lt = Lt;
    f2.norm_L = 4*sum((vector_tau(1:degree).*max(reshape(matrix_of_weights,N-1,degree))).^2); % square or not!? -> square

    start_point=sol2;
end
time_end_CA_hh=toc;

info2.objective=obj_f;
%% save into batch
batch{batch_cnt}.algorithm = 'Condat_hh_step';
batch{batch_cnt}.fdata_smooth = f1;
batch{batch_cnt}.f_nonsmooth = f2;
batch{batch_cnt}.sol = sol2.primal;
batch{batch_cnt}.info = info2;
batch{batch_cnt}.y_recon = y_recon;
batch{batch_cnt}.final_objective = f1.eval(f1.L(sol2.primal)) + f2.eval(f2.L(sol2.primal));
batch{batch_cnt}.elapsed_time = time_end_CA_hh;
batch{batch_cnt}.weights = reshape(matrix_of_weights,N-1,degree);

