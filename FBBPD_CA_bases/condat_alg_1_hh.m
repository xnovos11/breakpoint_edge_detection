function [sol, info] = condat_alg_1_hh(start_point, h1, h2, param)
% Condat algorithm 

if isa(start_point,'numeric') %single input, generate default dual variables
    temp = start_point;
    clear start_point
    start_point.primal = temp;
    start_point.dual1 = h1.L(start_point.primal);
    start_point.dual2 = h2.L(start_point.primal);
end

sol.primal = start_point.primal;
sol.dual1 = start_point.dual1;
sol.dual2 = start_point.dual2;

%step-sizes
param.condat_sigma = 1 / sqrt(h1.norm_L^2 + h2.norm_L^2);
param.condat_tau = 1 / sqrt(h1.norm_L^2 + h2.norm_L^2);
param.condat_ro=1.99;

%% Algorithm
x_old = sol.primal;
u1_old = sol.dual1;
u2_old = sol.dual2;
info.objective=nan(param.maxit,1);
tStart=tic;
for i=1:param.maxit
    x_half=x_old-param.condat_tau*(h1.Lt(u1_old)+h2.Lt(u2_old));
   
    sol.primal=param.condat_ro*x_half+(1-param.condat_ro)*x_old;
    
    u1_half=prox_adjoint(u1_old+param.condat_sigma*h1.L(2*x_half-x_old), param.condat_sigma, h1); 
    
    sol.dual1=param.condat_ro*u1_half+(1-param.condat_ro)*u1_old;
    
    u2_half=prox_adjoint(u2_old+param.condat_sigma*h2.L(2*x_half-x_old), param.condat_sigma, h2);
    
    sol.dual2=param.condat_ro*u2_half+(1-param.condat_ro)*u2_old;
   
    x_old=sol.primal;
    u1_old=sol.dual1;
    u2_old=sol.dual2;

    info.objective(i)=h1.eval(h1.L(sol.primal))+h2.eval(h2.L(sol.primal));
end
tElapsed=toc(tStart);
info.iter=i;
info.time=tElapsed;
