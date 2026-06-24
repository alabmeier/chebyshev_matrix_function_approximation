function w = cheb_matvec(A, v, a, b, c)

% CHEB_MATVEC Approximate f(A)v using a Chebyshev expansion.
%
% Computes an approximation of the matrix-vector product f(A)v
% without explicitly forming f(A). The approximation is evaluated
% using a Chebyshev polynomial expansion and the three-term
% recurrence relation.
%
% Inputs:
%   A - Square matrix
%   v - Input vector
%   a - Lower bound on the spectrum of A
%   b - Upper bound on the spectrum of A
%   c - Vector of Chebyshev coefficients
%
% Output:
%   w - Approximation of f(A)v
%
% Notes:
%   The matrix A is first mapped from [a,b] to [-1,1].
%   The algorithm requires only matrix-vector multiplications,
%   making it efficient for large sparse matrices.
%
% Example:
%   w = cheb_matvec(A,v,0,4,c);

alpha = 2/(b-a);
beta = -(a+b)/(b-a);

A_tilde = alpha*A + beta*speye(size(A));

v_prev = v;
v_curr = A_tilde*v;

w = c(1)*v_prev + c(2)*v_curr;

for k = 2:length(c)-1

    v_next = 2*(A_tilde*v_curr) - v_prev;

    w = w + c(k+1)*v_next;

    v_prev = v_curr;
    v_curr = v_next;

end

end
