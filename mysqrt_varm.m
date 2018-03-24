function X = mysqrt_varm(A)
%mysqrt_varm computes the matrix square root of a variance-covariance
%matrix.
%X = MYSQRTM_VARM(A) is the principal square root of the matrix A, i.e. X*X = A.
%          
%X is the unique square root for which every eigenvalue has a nonnegative
%real part.  
%The eigen decomposition always finds a diagonal matrix for which to
%compute the square root because a variance covariance matrix is always
%symmetric and positive-definite so the eigen decomposition always exists.
%In other words, Schur decomposition is not needed.
%For computation of square root of a general matrix, see sqrtm.m

[Q,T]=eig(A);          % T contains the eigenvalues of A
                       % Q contains the eigenvectors of A

R = diag(sqrt(diag(T)));     % Square root always exists. 

X = Q*R*Q';

