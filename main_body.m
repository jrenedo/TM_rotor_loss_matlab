% This code shows the calculation of rotor losses in high speed PM machines
% author: Jaime Renedo Anglada (renedo.jaime@ieee.org)

im=1i; 
% Define geometry:

mu_0=4*pi*10.0^-7; % [m kg s^-2 A^-2]
p=2; % number of pole pairs

B_r_ext=2.553*10.0^-3; % external field [T]
space_order=9; % space order of the harmonic
q=p*space_order; % harmonic order



n_mech=65000; % mechanical speed [rpm]
% n_mech=1500; % mechanical speed [rpm]
k_time=1*6; % time order
f_1=n_mech/60*p; % [Hz]


N_div=1000.0; % number of divisions simplified method


mu_r=1.07;

% Conductivity:
sigma_1=6.7*10.0^6; % [S/m] 
sigma_2=0.77*10.0^6; % [S/m]
sigma_3=2.2*10.0^4; % [S/m] 
sigma_4=3*10.0^-15; % [S/m] 
sigma_5=3*10.0^-15; % [S/m] 


% test:
sigma_3=3*10.0^-15; % [S/m]
%sigma_3=8*10^5
%sigma_sleeve=3*10.0^-15; % [S/m] 

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

r_wave=R_2+0.25*10.0^-3; % where you measured the harmonics

%B_r_ext=B_r_ext*(R_3/(R_2+0.5*10.0^-3))^(q-1);


L=109*10.0^-3; % [m] axial length


% no load losses

% vec_Br=[6.073 9.436 2.553 1.42 2.744 4.267 1.155]*10.0^-3; % vector with the amplitude of the asynchronous harmonics
vec_Br=[5.4033 8.703 2.4408 1.407 2.1913 3.5324 0.9911]*10.0^-3; % vector with the amplitude of the asynchronous harmonics
vec_so=[5 7 9 11 11 13 15]; % vector with the space orders of the asynchronous harmonics
vec_time=[1 1 1 1 2 2 2]; % vector with the time orders of the asynchronous harmonics


% Loop to calculate the power loss of each harmonic
vec_loss=zeros(1,length(vec_Br));
for count=1:length(vec_Br)
    space_order=vec_so(count)
    B_r_ext=vec_Br(count)
    k_time=vec_time(count)
    
    vec_loss(count)=Calc_rotor_loss( n_mech, space_order, p, k_time, B_r_ext,L,R_1,R_2,R_3,R_4,R_5,sigma_1,sigma_2,sigma_3,sigma_4,sigma_5,mu_r,r_wave);

    
end

count=1;

vec_loss_nl=vec_loss;



% on load losses

% delta=2/3
vec_Br=[6.073+6.79 9.436+4.21 2.553 1.42 2.744+1.622 4.267+1 1.155]*10.0^-3; % vector with the amplitude of the asynchronous harmonics
% delta=1
vec_Br=[5.4033 8.703 2.4408 1.407 2.1913 3.5324 0.9911]*10.0^-3; % vector with the amplitude of the asynchronous harmonics
vec_Br=[5.4033+7.8406 8.703+4.86 2.553 1.407 2.1913+1.827 3.5324+1 0.9911]*10.0^-3; % vector with the amplitude of the asynchronous harmonics

% old:
% vec_Br=[6.073+7.8406 9.436+4.86 2.553 1.42 2.744+1.827 4.267+1 1.1626]*10.0^-3; % vector with the amplitude of the asynchronous harmonics

vec_so=[5 7 9 11 11 13 15]; % vector with the space orders of the asynchronous harmonics
vec_time=[1 1 1 1 2 2 2]; % vector with the time orders of the asynchronous harmonics



% Loop to calculate the power loss of each harmonic
vec_loss=zeros(1,length(vec_Br));
for count=1:length(vec_Br)
    space_order=vec_so(count)
    B_r_ext=vec_Br(count)
    k_time=vec_time(count)
    
    vec_loss(count)=Calc_rotor_loss( n_mech, space_order, p, k_time, B_r_ext,L,R_1,R_2,R_3,R_4,R_5,sigma_1,sigma_2,sigma_3,sigma_4,sigma_5,mu_r,r_wave);

    
end

count=1;

vec_loss_ol=vec_loss;

vec_loss_nl
tot_loss_nl=sum(vec_loss_nl)
vec_loss_ol
tot_loss_ol=sum(vec_loss_ol)


