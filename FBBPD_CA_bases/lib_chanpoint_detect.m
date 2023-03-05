function  [breakpoints, diff_parameterization_coeff_filtred, diff_parameterization_coeff]=lib_chanpoint_detect(sol2, N, degree, param_method) 
 
 % Function that detects borders of segments in given parameterization
 % coefficients.
 %
 % INPUT:
 % sol2 - parameterization coefficients
 % N - length of signal y
 % degree - polynom degree
 % param_method - structure of paramters
 % - threshold - threshold used for breakpoint detection
 % - window_length - window length of median filter
 % - p - type of norm. Default p=2 -> Euclidean norm. 
 % - proximity - smallest proximity of two breakpoints. Default proximity=2
 %
 % OUTPUT:
 % breakpoints - detected breakpoints
 % diff_parameterization_coeff_filtred - vector of differences in
 % filtered Euclidean distance in parameterization coefficients
 %-------------------------------------------------------------------------

 threshold = param_method.threshold;
 window_length = param_method.window_length;
 parameterization_coeff_sol = reshape(sol2,N,degree); % reshape given parameters
 matrix_of_differences = zeros(N-1,degree);
 matrix_of_weighted_differences = zeros(N-1,degree);

% differences in paramteres
  for i = 1:degree
      matrix_of_differences(:,i) = diff(parameterization_coeff_sol(:,i));
  end

% Euclidean distance in differneces of parameters
 for i = 1:degree
     max_abs=(max(abs(matrix_of_differences(:,i))));
     if max_abs>0
        matrix_of_weighted_differences(:,i) = abs((matrix_of_differences(:,i))/max_abs);
     end
  %   disp(['alpha_', num2str(i-1) , ': ', num2str(1/(max(abs(matrix_of_differences(:,i)))))])
 end
 
% matrix_of_weighted_differences=matrix_of_differences;
 if ~isfield(param_method, 'p'), param_method.p = 2 ; end
 p = param_method.p;
 diff_parameterization_coeff = (sum(matrix_of_weighted_differences.^p,2)).^(1/p);

% Median filter of Euclidean distance of differences of parameters
median_diff = medfilt1(diff_parameterization_coeff,window_length);
diff_parameterization_coeff_filtred = diff_parameterization_coeff-median_diff;

% Breakpoint detection 
breakpoints = find(abs(diff_parameterization_coeff_filtred)>threshold); % find breakpoints, when the differences are larger than threshold
if length(breakpoints)>0    
    if breakpoints(1)==1
        breakpoints=breakpoints(2:end);
    end
    if ~isempty(breakpoints) && breakpoints(end)==N-1
        breakpoints=breakpoints(1:end-1);
    end
end
breakpoints = [0;breakpoints;N]; % defines start position of first segment minus 1, and end position of last segment

% delete breakpoint which are next to each other, leaves the one with
% larger value in difference vector
if ~isfield(param_method, 'proximity'), param_method.proximity = 2 ; end
proximity = param_method.proximity;
diff_breakpoints = (diff(breakpoints)<proximity);
pozice = find(diff_breakpoints==1);    
while sum(diff_breakpoints)>0
    for j = 1:length(pozice)
        i = pozice(j);
        if diff_breakpoints(i)==1 && breakpoints(i)<(N-(proximity-1))
            br_1 = diff_parameterization_coeff_filtred(breakpoints(i)); 
            br_2 = diff_parameterization_coeff_filtred(breakpoints(i+1));
            if abs(br_1)>abs(br_2)
                breakpoints_pom = [breakpoints(1:i);breakpoints(i+2:end)];
            else
                breakpoints_pom = [breakpoints(1:i-1);breakpoints(i+1:end)];
            end
            break
        elseif breakpoints(i)==0
            breakpoints_pom = [breakpoints(1:i);breakpoints(i+2:end)];
        end

    end
    breakpoints = breakpoints_pom;
    diff_breakpoints = (diff(breakpoints)<proximity);
    pozice = find(diff_breakpoints==1); 
end

%% plot of subresults

figure('Name','Differences in parameters..','NumberTitle','off')
 for i=1:degree
     matrix_of_differences(:,i)=diff(parameterization_coeff_sol(:,i));
     subplot (degree+2,1,i) 
     plot(matrix_of_differences(:,i))
     title (['Differences in parameter x_', num2str(i-1)])
 end
 
 subplot(degree+2,1,degree+1)
 plot(diff_parameterization_coeff)
 title ('Euclidean distance of differences in parameters (d)')
 
 % Median filter of Euclidean distance of differences of parameters
subplot (degree+2,1,degree+2) 
median_diff=medfilt1(diff_parameterization_coeff,window_length);
diff_parameterization_coeff_filtred=diff_parameterization_coeff-median_diff;
plot (diff_parameterization_coeff_filtred)
title (['Subtracted median filtering of d (with window of length ', num2str(window_length), ')'])