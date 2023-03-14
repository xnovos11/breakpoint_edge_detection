function [sol, info] = condat_alg_1_hhh(start_point, h1, h2, h3, param)
% Condat algorithm 1

if isa(start_point,'numeric') %single input, generate default dual variables
    temp = start_point;
    clear start_point
    start_point.primal = temp;
    start_point.dual1 = h1.L(start_point.primal);
    start_point.dual2 = h2.L(start_point.primal);
    start_point.dual3 = h3.L(start_point.primal);
end

sol.primal = start_point.primal;
sol.dual1 = start_point.dual1;
sol.dual2 = start_point.dual2;
sol.dual3 = start_point.dual3;

%step-sizes
param.condat_sigma = 1 / sqrt(h1.norm_L + h2.norm_L + h3.norm_L);
param.condat_tau = 1 / sqrt(h1.norm_L + h2.norm_L + h3.norm_L);

param.condat_ro=1.99;

%% Algorithm
x_old = sol.primal;
u1_old = sol.dual1;
u2_old = sol.dual2;
u3_old = sol.dual3;

info.objective=nan(param.maxit,1);

tStart=tic;

for i=1:param.maxit
    x_half=x_old-param.condat_tau*(h1.Lt(u1_old)+h2.Lt(u2_old)+h3.Lt(u3_old));
    
    sol.primal=param.condat_ro*x_half+(1-param.condat_ro)*x_old;
    
    p1 = u1_old+param.condat_sigma*h1.L(2*x_half-x_old);
    
    u1_half = p1 - param.condat_sigma*h1.prox(p1/param.condat_sigma, 1/param.condat_sigma);
       
    sol.dual1=param.condat_ro*u1_half+(1-param.condat_ro)*u1_old;
    
    p2 = u2_old+param.condat_sigma*h2.L(2*x_half-x_old);
    
    u2_half = p2 - param.condat_sigma*h2.prox_nov(p2/param.condat_sigma, 1/param.condat_sigma);
        
    sol.dual2=param.condat_ro*u2_half+(1-param.condat_ro)*u2_old;
    
    p3 = u3_old+param.condat_sigma*h3.L(2*x_half-x_old);
    
    u3_half = p3 - param.condat_sigma*h3.prox_nov(p3/param.condat_sigma, 1/param.condat_sigma);
    
    sol.dual3=param.condat_ro*u3_half+(1-param.condat_ro)*u3_old;
   
    x_old=sol.primal;
    u1_old=sol.dual1;
    u2_old=sol.dual2;
    u3_old=sol.dual3;

    info.objective(i)=h1.eval(h1.L(sol.primal))+h2.eval(h2.L(sol.primal))+h3.eval(h3.L(sol.primal));
    
end

tElapsed=toc(tStart);
info.iter=i;
info.time=tElapsed;
