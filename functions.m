% This code shows the functions to calculate the transfer matrices in high speed PM machines
% author: Jaime Renedo Anglada (renedo.jaime@ieee.org)


function TransferMatrix_comp( R_1, R_2, mu, sigma, omega, q)
% Calculation of the transfer matrix of a region. Complete methodology.
%   Detailed explanation goes here

    mu_0=4.0*pi*10.0^-7; % [m kg s^-2 A^-2]
    k_p=sqrt(im*omega*mu*sigma);

    F_2=mu_0*im*k_p*q/(mu*R_1)*(besselk(q,k_p*R_1)*besseli_d(q,k_p*R_1)-besseli(q,k_p*R_1)*besselk_d(q,k_p*R_1));
    M_2=[-mu_0*k_p/mu*besselk_d(q,k_p*R_1) -im*q/R_1*besselk(q,k_p*R_1); mu_0*k_p/mu*besseli_d(q,k_p*R_1) im*q/R_1*besseli(q,k_p*R_1)];
    N_2=[im*q/R_2*besseli(q,k_p*R_2) im*q/R_2*besselk(q,k_p*R_2); -mu_0*k_p/mu*besseli_d(q,k_p*R_2) -mu_0*k_p/mu*besselk_d(q,k_p*R_2)];
    T_mat=[1.0 0.0; 0.0 -1.0]*N_2*M_2/F_2*[1.0 0.0; 0.0 -1.0];
    
return T_mat
end

function besseli_d(q,x)
% Derivative of the modified Bessel function.
%   Detailed explanation goes here

    I_prima=besseli(q-1,x)-q/x*besseli(q,x)
    
return I_prima
end

function besselk_d(q,x)
% Derivative of the modified Bessel function.
%   Detailed explanation goes here

    K_prima=-besselk(q-1,x)-q/x*besselk(q,x)
    
return K_prima
end

function TransferMatrix_medium_simp( R_1, R_2, mu, sigma, omega, q, N_div)
% Calculation of the transfer matrix of a region. Medium simplified methodology.
%   Detailed explanation goes here

mu_0=4.0*pi*10.0^-7; % [m kg s^-2 A^-2]

S_2l=R_2-R_1;
S_2d=S_2l/N_div;

T_2=[1.0 0.0; 0.0 1.0];
for count2 in 0:(N_div)
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
return T_2
end

function TransferMatrix_simp( R_1, R_2, mu, sigma, omega, q, N_div)
% Calculation of the transfer matrix of a region. Medium simplified methodology.
%   Detailed explanation goes here

mu_0=4.0*pi*10.0^-7; % [m kg s^-2 A^-2]

S_2l=R_2-R_1;
S_2d=S_2l/N_div;

T_2=[1.0 0.0; 0.0 1.0];
for count2 in 0:(N_div)
    delta=S_2d;
    r_2a=R_1+delta*count2+delta/2;
    k_2=q/r_2a;
    d_2=1/sqrt(sigma*mu*omega);
    gamma_2h=sqrt(k_2^2+im/d_2^2);
    beta_2h=gamma_2h/(im*mu*k_2);
    T_2h=[1.0 gamma_2h*S_2d/beta_2h/mu_0; mu_0*beta_2h*gamma_2h*S_2d 1.0];
%     T_2h=[cosh(gamma_2h*S_2d) sinh(gamma_2h*S_2d)/beta_2h/mu_0; mu_0*beta_2h*sinh(gamma_2h*S_2d) cosh(gamma_2h*S_2d)];
    
    T_2=T_2*T_2h;
    
end
return T_2
end

function TransferMatrix_super_simp( R_1, R_2, mu, sigma, omega, q, N_div)
%Calculation of the transfer matrix of a region. Medium simplified methodology.
%   Detailed explanation goes here

mu_0=4.0*pi*10.0^-7; % [m kg s^-2 A^-2]

S_2l=R_2-R_1;
S_2d=S_2l/N_div;
    
delta=S_2d;
r_2a=(R_1+R_2)/2;
k_2=q/r_2a;
d_2=1/sqrt(sigma*mu*omega);
gamma_2h=sqrt(k_2^2+im/d_2^2);
beta_2h=gamma_2h/(im*mu*k_2);
    
T_2h=[1.0 gamma_2h*S_2d/beta_2h/mu_0; mu_0*beta_2h*gamma_2h*S_2d 1.0];
%     T_2h=[cosh(gamma_2h*S_2d) sinh(gamma_2h*S_2d)/beta_2h/mu_0; mu_0*beta_2h*sinh(gamma_2h*S_2d) cosh(gamma_2h*S_2d)];

T_2=T_2h^N_div;
    
return T_2
end

function TransferMatrix_nc( R_1, R_2, mu, sigma, omega, q)
% Calculation of the transfer matrix of a region. Complete methodology.
%   Detailed explanation goes here

    mu_0=4.0*pi*10.0^-7; % [m kg s^-2 A^-2]
    k_p=sqrt(im*omega*mu*sigma);

    F_2=mu_0*im*q/(mu*R_1*R_2);
    M_2=[-mu_0*R_1^(-q)/mu -im*q*R_1^(-q-1); mu_0*R_1^(q)/mu im*q*R_1^(q-1)];
    N_2=[im*q*R_1^(q-1) im*q*R_1^(-q-1); -mu_0*R_1^(q)/mu -mu_0*R_1^(-q)/mu];
    T_mat=[1.0 0.0; 0.0 -1.0]*N_2*M_2/F_2*[1.0 0.0; 0.0 -1.0];
    
return T_mat
end



