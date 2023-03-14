function Y = soft_thresh_group(X,T)

% X is a matrix where rows are the groups
% Y is the group-thresholded output of the same size
% T is the threshold value

% number of columns
N = size(X,2);

% compute norms of the rows
norms = sqrt(sum(abs(X).^2,2));

norms(norms==0)=eps;

% soft threshold
Y = X./repmat(norms,1,N).* max(repmat(norms-T,1,N),0);