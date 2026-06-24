function c = cheb_coeff(f, a, b, n, ng)

% CHEB_COEFF Compute Chebyshev expansion coefficients.
%
% Computes the first n+1 Chebyshev coefficients for approximating
% a function f on the interval [a,b] using Gauss-Chebyshev
% quadrature with ng nodes.
%
% Inputs:
%   f  - Function handle for the target function
%   a  - Left endpoint of approximation interval
%   b  - Right endpoint of approximation interval
%   n  - Degree of Chebyshev approximation
%   ng - Number of Gauss-Chebyshev quadrature points
%
% Output:
%   c  - Vector of Chebyshev coefficients
%
% Notes:
%   The interval [a,b] is mapped to [-1,1] before computing
%   coefficients. Coefficients are obtained using numerical
%   integration with Gauss-Chebyshev quadrature.
%
% Example:
%   f = @(x) exp(-0.01*x);
%   c = cheb_coeff(f,0,100,50,100);

k = 1:ng;

x = cos((2*k-1)*pi/(2*ng));
w = pi/ng * ones(size(x));

x_ab = ((b-a)*x + (b+a))/2;

f_eval = f(x_ab);

c = zeros(n+1,1);

for j = 0:n

    Tj = cheb_eval(j,x);

    c(j+1) = (2/(pi*(1+(j==0)))) * ...
        sum(w .* f_eval .* Tj);

end

end
