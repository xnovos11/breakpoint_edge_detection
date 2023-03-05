function [sol,info] = proj_ball(x, param)
% Projection on the ball
% x
% param.y - data
% param.epsilon

% sol

%% 
sol=(param.epsilon*(x-param.y))/max(norm(x-param.y),param.epsilon)+param.y;
if sum(x-sol)==0
    info.change='No';
else
    info.change='Yes';
end