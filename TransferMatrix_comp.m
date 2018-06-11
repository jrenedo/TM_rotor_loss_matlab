function T_mat=TransferMatrix_comp( R_1, R_2, mu, sigma, omega, q)
% Calculation of the transfer matrix of a region. Complete methodology.
%   Detailed explanation goes here
    
    im=1i; 
    mu_0=4.0*pi*10.0^-7; % [m kg s^-2 A^-2]
    k_p=sqrt(im*omega*mu*sigma);

    F_2=mu_0*im*k_p*q/(mu*R_1)*(besselk(q,k_p*R_1)*besseli_d(q,k_p*R_1)-besseli(q,k_p*R_1)*besselk_d(q,k_p*R_1));
    M_2=[-mu_0*k_p/mu*besselk_d(q,k_p*R_1) -im*q/R_1*besselk(q,k_p*R_1); mu_0*k_p/mu*besseli_d(q,k_p*R_1) im*q/R_1*besseli(q,k_p*R_1)];
    N_2=[im*q/R_2*besseli(q,k_p*R_2) im*q/R_2*besselk(q,k_p*R_2); -mu_0*k_p/mu*besseli_d(q,k_p*R_2) -mu_0*k_p/mu*besselk_d(q,k_p*R_2)];
    T_mat=[1.0 0.0; 0.0 -1.0]*N_2*M_2/F_2*[1.0 0.0; 0.0 -1.0];
    

end