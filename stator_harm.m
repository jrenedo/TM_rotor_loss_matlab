% This code shows the calculation of rotor losses in high speed PM machines
% author: Jaime Renedo Anglada (renedo.jaime@ieee.org)


% Define geometry:

mu_0=4*pi*10.0^-7; % [m kg s^-2 A^-2]

p=2; % number of pole pairs
Qs=12.0; % number of slots
N_ph=3; % number of phases

N_turns=7; % number of turns per phase
I_n=60; % current [A]

m_spp=Qs/(2.0*p*N_ph); % number of slots per pole per phase
beta_pa=2*pi/Qs; % slot pitch angle

%delta_p=2/3; % ratio of short corded to full corded coil pitch
delta_p=1; % ratio of short corded to full corded coil pitch
% What value shall I put here? 1 right? you mentioned that it is fully pitched. 
alpha_s=0; % skew angle

k_space=1; % space order of the harmonic
q=p*k_space; % harmonic order



n_mech=65000; % mechanical speed [rpm]
k_time=1; % time order
f_1=n_mech/60*p; % electrical frequency [Hz]
omega=2*pi*f_1; % electrical frequency [rad/s]


mu_r=1.07;

% Conductivity:
sigma_1=6.7*10.0^6; % [S/m] 
sigma_2=0.77*10.0^6; % [S/m]
sigma_3=2.2*10.0^4; % [S/m] 
sigma_4=3*10.0^-15; % [S/m] 
sigma_5=3*10.0^-15; % [S/m] 



% Permeability:
mu_1=750*mu_0; % [m kg s^-2 A^-2]
mu_2=mu_r*mu_0; % [m kg s^-2 A^-2]
mu_3=1*mu_0; % [m kg s^-2 A^-2]
mu_4=1*mu_0; % [m kg s^-2 A^-2]
mu_5=5000*mu_0; % [m kg s^-2 A^-2]



% Geometry:
t_sl=2.95*10.0^-3; % [m] sleeve thickness

R_1=20*10.0^-3; % [m] rotor radius
R_2=27.4*10.0^-3; % [m] PM radius
R_3=R_2+t_sl; % [m] sleeve outer radius
R_4=31.15*10.0^-3; % [m] current sheet radius
R_5=50*10.0^-3; % [m] end domain radius

% where to evaluate the harmonics:
r_wave=R_2+0.25*10^-3; % where to calculate the harmonics
% r_wave=R_3; % where to calculate the harmonics

L=109*10.0^-3; % [m] axial length





% Calculation of the winding factors:

K_dq=sin(m_spp*beta_pa*p*k_space/2)/(m_spp*sin(beta_pa*p*k_space/2));
K_dq

K_pq=sin(k_space*delta_p*pi/2);
K_pq

%K_sq=sin(k_space*alpha_s/2)/(k_space*alpha_s/2);
K_sq=1;
K_sq

K_wq=K_dq*K_pq*K_sq;
K_wq



% Calculation of the MMF:
k_space


F_qk=3/2*4/pi*N_turns/(2*p)*1/q*K_wq*2.0^(0.5)*I_n; % [A]
%F_qk=F_qk*2


% Calculation of the amplitude of the harmonics:

J_qk=q*p/R_4*F_qk; % in our case R_4 is the stator bore

B_rq=J_qk*R_4*mu_0*(R_1^(-2*q)*r_wave^q+r_wave^-q)/(r_wave*(R_1^(-2*q)*R_4^q-R_4^-q)); % [T]
B_tq=J_qk*R_4*mu_0*(R_1^(-2*q)*r_wave^q-r_wave^-q)/(r_wave*(R_1^(-2*q)*R_4^q-R_4^-q)); % [T]

B_rq*1000
B_tq*1000


