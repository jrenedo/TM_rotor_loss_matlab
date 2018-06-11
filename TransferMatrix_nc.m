function T_mat=TransferMatrix_nc( R_1, R_2, mu, sigma, omega, q)
% Calculation of the transfer matrix of a region. Complete methodology.
%   Detailed explanation goes here
    
    im=1i; 
    mu_0=4.0*pi*10.0^-7; % [m kg s^-2 A^-2]
    k_p=sqrt(im*omega*mu*sigma);

    F_2=mu_0*im*q/(mu*R_1*R_2);
    M_2=[-mu_0*R_1^(-q)/mu -im*q*R_1^(-q-1); mu_0*R_1^(q)/mu im*q*R_1^(q-1)];
    N_2=[im*q*R_1^(q-1) im*q*R_1^(-q-1); -mu_0*R_1^(q)/mu -mu_0*R_1^(-q)/mu];
    T_mat=[1.0 0.0; 0.0 -1.0]*N_2*M_2/F_2*[1.0 0.0; 0.0 -1.0];
    

end