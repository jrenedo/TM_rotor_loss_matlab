function [ I_prima ] = besseli_d(q,x)
% Derivative of the modified Bessel function.
%   Detailed explanation goes here
I_prima=besseli(q-1,x)-q/x*besseli(q,x);

end

