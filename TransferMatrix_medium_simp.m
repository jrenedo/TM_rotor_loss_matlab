function T_2=TransferMatrix_medium_simp( R_1, R_2, mu, sigma, omega, q, N_div)
% Calculation of the transfer matrix of a region. Medium simplified methodology.
%   Detailed explanation goes here
    
    im=1i; 
    mu_0=4.0*pi*10.0^-7; % [m kg s^-2 A^-2]
    S_2l=R_2-R_1;
S_2d=S_2l/N_div;

T_2=[1.0 0.0; 0.0 1.0];
for count2 = 0:(N_div)
    delta=S_2d;
    r_2a=R_1+delta*count2+delta/2;
    k_2=q/r_2a;
    d_2=1/sqrt(sigma*mu*omega);
    gamma_2h=sqrt(k_2^2+im/d_2^2);
    beta_2h=gamma_2h/(im*mu*k_2);
%    T_2h=[1.0 gamma_2h*S_2d/beta_2h/mu_0; mu_0*beta_2h*gamma_2h*S_2d 1.0];
    T_2h=[cosh(gamma_2h*S_2d) sinh(gamma_2h*S_2d)/beta_2h/mu_0; mu_0*beta_2h*sinh(gamma_2h*S_2d) cosh(gamma_2h*S_2d)];
    
    T_2=T_2*T_2h;
    
end

end

