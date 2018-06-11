% This code shows the functions to calculate the rotor losses in high speed PM machines
% author: Jaime Renedo Anglada (renedo.jaime@ieee.org)


function Calc_rotor_loss( n_mech, space_order, p, k_time, B_r_ext,L,R_1,R_2,R_3,R_4,R_5,sigma_1,sigma_2,sigma_3,sigma_4,sigma_5,mu_r,r_wave)
% Calculation of the rotor losses in one function

    % define geometry:
    mu_0=4*pi*10.0^-7; % [m kg s^-2 A^-2]
    q=p*space_order; % harmonic order
    k_time=1*6; % time order
    f_1=n_mech/60*p; % [Hz]
    f=f_1;


    % Permeability:
    mu_1=750*mu_0; % [m kg s^-2 A^-2]
    mu_2=mu_r*mu_0; % [m kg s^-2 A^-2]
    mu_3=1*mu_0; % [m kg s^-2 A^-2]
    mu_4=1*mu_0; % [m kg s^-2 A^-2]
    mu_5=5000*mu_0; % [m kg s^-2 A^-2]

    B_r_ext=B_r_ext*(R_3/r_wave)^(q-1); % where you measured the amplitude of the harmonics

    %%% No eddy currents:
    f_no_eddy=0.01;
    omega=2*pi*f_no_eddy; % [rad/s] frequency
    J_kq=1; % current sheet density
    sigma_air=3*10.0^-15; % [S/m] 

    % Propagation constants:
    k_p_1=sqrt(im*omega*mu_1*sigma_air);
    k_p_2=sqrt(im*omega*mu_2*sigma_air);
    k_p_3=sqrt(im*omega*mu_3*sigma_air);
    k_p_4=sqrt(im*omega*mu_4*sigma_air);
    k_p_5=sqrt(im*omega*mu_5*sigma_air);


    % Region 1:
    %beta_1=mu_0*im*k_p_1*R_1*besseli_d(q,k_p_1*R_1)/besseli(q,k_p_1*R_1)/mu_1/q;
    beta_1=mu_0*im/mu_1/q;


    % Region 2:
    %T_2=TransferMatrix_nc( R_1, R_2, mu_2, sigma_air, omega, q);
    T_2=TransferMatrix_simp( R_1, R_2, mu_2, sigma_air, omega, q, N_div);
    %         T_2_pm=TransferMatrix_medium_simp( R_1, R_2, mu_2, sigma_2, omega, q, N_div);
    %         T_2_pm=TransferMatrix_simp( R_1, R_2, mu_2, sigma_2, omega, q, N_div);
    %         T_2_pm=TransferMatrix_super_simp( R_1, R_2, mu_2, sigma_2, omega, q, N_div);
  

    % Region 3:
    %T_3=TransferMatrix_nc( R_2, R_3, mu_3, sigma_air, omega, q);
    T_3=TransferMatrix_simp( R_2, R_3, mu_3, sigma_air, omega, q, N_div);


    % Region 4:
    %T_4=TransferMatrix_nc( R_3, R_4, mu_4, sigma_air, omega, q);
    T_4=TransferMatrix_simp( R_3, R_4, mu_4, sigma_air, omega, q, N_div);


    % Region 5:
    %T_5=TransferMatrix_nc( R_4, R_5, mu_5, sigma_air, omega, q);
    T_5=TransferMatrix_simp( R_4, R_5, mu_5, sigma_air, omega, q, N_div);



    % Estimation of the coefficients:

    Mat_D=[0.0 0.0; 0.0 1.0]-T_5*T_4*T_3*T_2*[1.0 0.0; beta_1 0.0];

    %         temp_mat=inv(Mat_D);
    vec_sol=inv(Mat_D)*T_5*[0; mu_0*J_kq];

    % Get the fields:
    B_1=vec_sol(1);
    H_1=-beta_1*B_1/mu_0;



    H_5=vec_sol(2)/mu_0;

    temp_2=T_2*[B_1; mu_0*(H_1)];
    B_2=temp_2(1);
    H_2=temp_2(2)/mu_0;



    temp_3=T_3*[B_2; mu_0*(H_2)];
    B_3=temp_3(1);
    H_3=temp_3(2)/mu_0;
  
  

    temp_4=T_4*[B_3; mu_0*(H_3)];
    B_4=temp_4(1);
    H_4=temp_4(2)/mu_0;


    % Field at the rotor hub:
    vec_hub_1=real(B_1+mu_1*H_1);
    vec_hub_1=abs(B_1);

    vec_hub_2=imag(B_1+mu_1*H_1);
    vec_hub_2=abs(mu_1*H_1);

    B_r_hub_ne=vec_hub_1;
    B_theta_hub_ne=vec_hub_2;

    % Field at the outer part of the PM:
    vec_out_1=real(B_2+mu_2*H_2);  
    vec_out_1=abs(B_2);

    vec_out_2=imag(B_2+mu_2*H_2);
    vec_out_2=abs(mu_2*H_2);

    B_r_pm_ne=vec_out_1;
    B_theta_pm_ne=vec_out_2;

    % Field at the outer part of the sleeve:
    vec_out_1=real(B_3+mu_3*H_3);  
    vec_out_1=abs(B_3);

    vec_out_2=imag(B_3+mu_3*H_3);
    vec_out_2=abs(mu_3*H_3);

    B_r_sleeve_ne=vec_out_1;
    B_theta_sleeve_ne=vec_out_2;

    
    %%% with eddy currents:
    f_2=f
    q=q;

    omega=2*pi*f_2*k_time; % [rad/s] frequency



    % Propagation constants:
    k_p_1=sqrt(im*omega*mu_1*sigma_1);
    k_p_2=sqrt(im*omega*mu_2*sigma_2);
    k_p_3=sqrt(im*omega*mu_3*sigma_3);
    k_p_4=sqrt(im*omega*mu_4*sigma_4);
    k_p_5=sqrt(im*omega*mu_5*sigma_5);


    % Region 1:
    beta_1=mu_0*im*k_p_1*R_1*besseli_d(q,k_p_1*R_1)/besseli(q,k_p_1*R_1)/mu_1/q;
    beta_1=mu_0*im/mu_1/q;

    % Region 2:
    T_2=TransferMatrix_comp( R_1, R_2, mu_2, sigma_2, omega, q);
    %         T_2_pm=TransferMatrix_medium_simp( R_1, R_2, mu_2, sigma_2, omega, q, N_div);
    %         T_2_pm=TransferMatrix_simp( R_1, R_2, mu_2, sigma_2, omega, q, N_div);
    %         T_2_pm=TransferMatrix_super_simp( R_1, R_2, mu_2, sigma_2, omega, q, N_div);

    % Region 3:
    T_3=TransferMatrix_comp( R_2, R_3, mu_3, sigma_3, omega, q);

    % Region 4:
    T_4=TransferMatrix_comp( R_3, R_4, mu_4, sigma_4, omega, q);


    % Region 5:
    T_5=TransferMatrix_comp( R_4, R_5, mu_5, sigma_5, omega, q);


    % Estimation of the coefficients:

    Mat_D=[0.0 0.0; 0.0 1.0]-T_5*T_4*T_3*T_2*[1.0 0.0; beta_1 0.0];
    %         temp_mat=inv(Mat_D);
    vec_sol=inv(Mat_D)*T_5*[0; mu_0*J_kq];
    % Get the fields:
    B_1=vec_sol(1);
    H_1=-beta_1*B_1/mu_0;


    H_5=vec_sol(2)/mu_0;


    temp_2=T_2*[B_1; mu_0*(H_1)];
    B_2=temp_2(1);
    H_2=temp_2(2)/mu_0;


    temp_3=T_3*[B_2; mu_0*(H_2)];
    B_3=temp_3(1);
    H_3=temp_3(2)/mu_0;


    temp_4=T_4*[B_3; mu_0*(H_3)];
    B_4=temp_4(1);
    H_4=temp_4(2)/mu_0;


    % Field at the rotor hub:
    vec_hub_1=real(B_1+mu_1*H_1);
    vec_hub_1=abs(B_1);

    vec_hub_2=imag(B_1+mu_1*H_1);
    vec_hub_2=abs(mu_1*H_1);

    B_r_hub=vec_hub_1;
    B_theta_hub=vec_hub_2;

    % Field at the outer part of the PM:
    vec_out_1=real(B_2+mu_2*H_2);  
    vec_out_1=abs(B_2);

    vec_out_2=imag(B_2+mu_2*H_2);
    vec_out_2=abs(mu_2*H_2);

    B_r_pm=vec_out_1;
    B_theta_pm=vec_out_2;

    % Field at the outer part of the sleeve:
    vec_out_1=real(B_3+mu_3*H_3);  
    vec_out_1=abs(B_3);

    vec_out_2=imag(B_3+mu_3*H_3);
    vec_out_2=abs(mu_3*H_3);

    B_r_sleeve=vec_out_1;
    B_theta_sleeve=vec_out_2;

    %%% Calculate the coefficients and power losses:

    S_1=2*pi*R_1*L; % area surf rotor steel
    S_2=2*pi*R_2*L; % area surf PM
    S_3=2*pi*R_3*L; % area surf sleeve
    S_4=2*pi*R_4*L; % area surf bore


    % region 2 (PM):
    mu=mu_2
    sigma=sigma_2
    q=q
    k_p=k_p_2;

    F_2=mu_0*im*k_p*q/(mu*R_1)*(besselk(q,k_p*R_1)*besseli_d(q,k_p*R_1)-besseli(q,k_p*R_1)*besselk_d(q,k_p*R_1));
    M_2=[-mu_0*k_p/mu*besselk_d(q,k_p*R_1) -im*q/R_1*besselk(q,k_p*R_1); mu_0*k_p/mu*besseli_d(q,k_p*R_1) im*q/R_1*besseli(q,k_p*R_1)];
    N_2=[im*q/R_2*besseli(q,k_p*R_2) im*q/R_2*besselk(q,k_p*R_2); -mu_0*k_p/mu*besseli_d(q,k_p*R_2) -mu_0*k_p/mu*besselk_d(q,k_p*R_2)];
    T_mat=[1.0 0.0; 0.0 -1.0]*N_2*M_2/F_2*[1.0 0.0; 0.0 -1.0];


    temp=M_2*[B_1; mu_0*H_1]/F_2;

    C_2=temp(1);
    D_2=temp(2);


    E_2=-im*omega*(C_2*besseli(q,k_p_2*R_2)+D_2*besselk(q,k_p_2*R_2)); % electric field intensity
    H_theta_2=-k_p_2*(C_2*besseli_d(q,k_p_2*R_2)+D_2*besselk_d(q,k_p_2*R_2))/mu_0; % H_theta complex value

    E_2=omega*R_2/q*B_2; % electric field intensity
    H_theta_2=H_2; % H_theta complex value




    % region 3 (sleeve):
    mu=mu_3
    sigma=sigma_3
    q=q
    k_p=k_p_3;

    F_3=mu_0*im*k_p*q/(mu*R_2)*(besselk(q,k_p*R_2)*besseli_d(q,k_p*R_2)-besseli(q,k_p*R_2)*besselk_d(q,k_p*R_2));
    M_3=[-mu_0*k_p/mu*besselk_d(q,k_p*R_2) -im*q/R_2*besselk(q,k_p*R_2); mu_0*k_p/mu*besseli_d(q,k_p*R_2) im*q/R_2*besseli(q,k_p*R_2)];
    N_3=[im*q/R_3*besseli(q,k_p*R_3) im*q/R_3*besselk(q,k_p*R_3); -mu_0*k_p/mu*besseli_d(q,k_p*R_3) -mu_0*k_p/mu*besselk_d(q,k_p*R_3)];
    T_mat=[1.0 0.0; 0.0 -1.0]*N_3*M_3/F_3*[1.0 0.0; 0.0 -1.0];


    temp=M_3*[B_2; mu_0*H_2]/F_3;

    C_3=temp(1);
    D_3=temp(2);


    E_3=-im*omega*(C_3*besseli(q,k_p_3*R_3)+D_3*besselk(q,k_p_3*R_3)); % electric field intensity
    H_theta_3=-k_p_3*(C_3*besseli_d(q,k_p_3*R_3)+D_2*besselk_d(q,k_p_3*R_3))/mu_0; % H_theta complex value

    E_3=omega*R_3/q*B_3; % electric field intensity
    H_theta_3=H_3; % H_theta complex value



    % Test of matrix:
    T_2_comp=T_3*T_2


    % region 4 (air-gap):
    mu=mu_4
    sigma=sigma_4
    q=q
    k_p=k_p_4;

    F_4=mu_0*im*k_p*q/(mu*R_3)*(besselk(q,k_p*R_3)*besseli_d(q,k_p*R_3)-besseli(q,k_p*R_3)*besselk_d(q,k_p*R_3));
    M_4=[-mu_0*k_p/mu*besselk_d(q,k_p*R_3) -im*q/R_3*besselk(q,k_p*R_3); mu_0*k_p/mu*besseli_d(q,k_p*R_3) im*q/R_3*besseli(q,k_p*R_3)];
    N_4=[im*q/R_4*besseli(q,k_p*R_4) im*q/R_4*besselk(q,k_p*R_4); -mu_0*k_p/mu*besseli_d(q,k_p*R_4) -mu_0*k_p/mu*besselk_d(q,k_p*R_4)];
    T_mat=[1.0 0.0; 0.0 -1.0]*N_3*M_3/F_3*[1.0 0.0; 0.0 -1.0];

    temp=M_4*[B_3; mu_0*H_3]/F_4;

    C_4=temp(1);
    D_4=temp(2);


    % Scale coefficient: 
    K_B=B_r_ext/B_r_sleeve_ne;


    % Power going towards the rotor sleeve:
    P3_temp=1/2*real(E_3*conj(H_theta_3))*S_3*(K_B)^2;
    Loss_factor_temp=P3_temp/B_r_ext^2/omega/S_3;

    % Power going towards the rotor PM:
    P2_temp=1/2*real(E_2*conj(H_theta_2))*S_2*(K_B)^2;

    % Power going towards the rotor hub:
    E_1=-im*omega*(C_2*besseli(q,k_p_2*R_1)+D_2*besselk(q,k_p_2*R_1)); % electric field intensity
    H_theta_1=-k_p_2*(C_2*besseli_d(q,k_p_2*R_1)+D_2*besselk_d(q,k_p_2*R_1))/mu_0; % H_theta complex value

    E_1=omega*R_1/q*B_1; % electric field intensity
    H_theta_1=H_1; % H_theta complex value

    P1_temp=1/2*real(E_1*conj(H_theta_1))*S_1*(K_B)^2;

    P_mag=P2_temp-P1_temp;
    P_hub=P1_temp;
    P_sleeve=P3_temp-P2_temp;
    Loss_factor=Loss_factor_temp;


return P3_temp
end