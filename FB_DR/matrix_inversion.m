function inversion = matrix_inversion(A,lambda) 

% Function which computes inversion of matrix I + lambda*A'A
% A = (I D D^2 D^3 ... D^k) or A = W(I D D^2 D^3 ... D^k)
% with scaling parameter W -- diagonal matrix of size D 
% A is recomended to be sparse!!!
% D is diagonal matrix (normalized to be max 1)
% ------------------------------------------------------------------------
[N,m] = size(A);
I = speye(N); 

c = I + lambda*(A*A');
inversion = speye(N*(m/N)) - lambda*kron(speye(m/N),sparse(1:N,1:N,1./diag(c)))*(A'*A);

end