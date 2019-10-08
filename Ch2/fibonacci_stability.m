function [L1, L2] = fibonacci_stability

% Illustrate (in)stability of Fibonacci's rabbit problem

% initialize parameters
T = 10;
N0 = 1;
N =   zeros(T,1);
N(1) = 0;
N(2) = N0;

% the solution to the characteristic equation, solved using quadratic
% formula:
lam1 = (1+5^0.5)/2;
lam2 = (1-5^0.5)/2;

% Given solution of dN(t) = c1*lam1^t + c2*lam2^t, solve for c1 and c2,
% given that we know dN(0) = 1 and dN(-1) = 0:
c2 = N0*(2-lam1)/( lam2* (lam2-lam1));
c1 = (N0 - c2 * lam2)/lam1;

for t = 2:(T+1)
    N(t+1) = N(t) + N(t-1);
    L1(t) = lam1^(t-2);
    L2(t) = lam2^(t-2);
end

