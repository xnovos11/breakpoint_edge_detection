function [batch, batch_cnt]=lib_signal_denoising(y, batch, batch_cnt, N, degree, param_method)
% Function that performs changepoints detection, least squres refitting of
% each detected segments and final signal reconstruction
%
% INPUT:
% batch - structure of results
% batch cnt - number of results in batch
% N - length of signal y
% degree - polynom degree
% param_method - structure of paramters
% - threshold - threshold used for breakpoint detection
% - window_length - window length of median filter
% - p - type of norm. Default p=2 -> Euclidean norm.
% - proximity - smallest proximity of two breakpoints. Default proximity=2
%
% OUTPUT:
% batch - structure of results
% batch cnt - number of results in batch
%--------------------------------------------------------------------------

for cnt = 1:batch_cnt
    % detection of breakpoints in given parameters
    [breakpoints, diff_param, diff_parameterization_coeff] = lib_chanpoint_detect(batch{cnt}.sol, N, degree, param_method);
    
    % saving the detected breakpoints in batch 
    batch{cnt}.breakpoints = breakpoints;
    batch{cnt}.diff_med_filt = diff_param;
    batch{cnt}.diff_euclidean = diff_parameterization_coeff;
    
    % refitting of each detected segment by using least squares method
    recon=zeros(N,degree);
    for i = 1:length(breakpoints)-1
        start = breakpoints(i)+1; %start position of segment
        stop = breakpoints(i+1); %stop position of segment
        
        % getting new paramtererization coefficients
        beta_det = lib_least_squares_fit(start, stop, y(start:stop), param_method);
        
        % creating new parameterizaion coefficients
        for j = 1:degree
            recon(start:stop,j) = ones(stop-start+1,1).*beta_det(j);
        end
    end
    
    % saving reconstruted signal and new parameters
    batch{cnt}.param_recon = recon;
    batch{cnt}.signal_recon = lib_recon_from_params(reshape(recon,degree*N,1),N,degree, param_method); % reconstructed signal
end

