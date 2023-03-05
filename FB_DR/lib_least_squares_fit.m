function solution = lib_least_squares_fit(start, stop, y, param_method)

% Function refitting a given segment of the signal by the ordinary least
% squares method, prosucining new parameters
%
% INPUT:
% start - segment start position in signal 
% stop - segment stop position in signal 
% y - observed signal - (column) vector of data
%
% OUTPUT: 
% solution ... M parameters minimizing square error ... [x0,x1, ... xM],
% number of parameres is given by degree
% -----------------------------------------------------------------------
matrix_of_d = param_method.A;

A = matrix_of_d([start:stop],:);
solution = A\y; % new parameters

