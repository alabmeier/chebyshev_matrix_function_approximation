function T = cheb_eval(n, x)

% CHEB_EVAL Evaluate the Chebyshev polynomial T_n(x).
%
% Computes the Chebyshev polynomial of the first kind of degree n
% at the points specified by x using the three-term recurrence
% relation.
%
% Inputs:
%   n - Nonnegative integer polynomial degree
%   x - Scalar, vector, or array of evaluation points
%
% Output:
%   T - Values of T_n(x) evaluated at x
%
% Example:
%   x = linspace(-1,1,100);
%   T3 = cheb_eval(3,x);

if n == 0
    T = ones(size(x));
elseif n == 1
    T = x;
else
    T_nminus1 = ones(size(x));
    T_n = x;

    for k = 2:n
        T_nplus1 = 2*x.*T_n - T_nminus1;
        T_nminus1 = T_n;
        T_n = T_nplus1;
    end

    T = T_n;
end

end
