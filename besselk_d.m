function [ K_prima ] = besselk_d(q,x)
% Derivative of the modified Bessel function.
%   Detailed explanation goes here
K_prima=-besselk(q-1,x)-q/x*besselk(q,x)

end