
MP_ON = 1;

% Max continuous
USE_MAX_OP = 0;
max_i = 10;
max_trq = 501e-3;

% Efficiencies
n_g = 1;

% gear ratio
N = 1;   

% acceleration due to gravity
g = 9.81;

angl = 0;
angl = deg2rad(angl);

% masses
m_p = 0.5;    % pulley  1  
m_b = 0.5;    % belt 1.5
m_s = 0.5;    % shaft 0.6
m_l = 25;    % load 25(no max)

% radii
r_p = 0.5e-1;    % pulley
r_s = 0.3e-1;    % shaft


% inertia
J_m = 0.000129;    % motor 
J_g = 0.0000;    % gearbox

J_p = m_p * r_p^2;  % pulley
J_l = m_l * r_p^2;  % load
J_b = m_b * r_p^2;  % belt
J_s = 0.5 * m_s * r_s^2;    % shaft

% Coefficients & Constants
b = 1e-5;      % shaft damping 1 * 10e-5
% mu_v = 87.38e-6;   % [viscous friction JUST FOR THE MOTOR, NOT THE LOAD!!! ADD LOAD] WRONG!!!  
b_m = 87.39e-6; % ----- THIS IS THE DAMPING EXPERIENCED BY THE MOTOR AT NO LOAD --- CORRECT

mu_s = 0.23;   % static friction (bed and belt) or couloumb?


Ke = 0.0537; % motor constant
R = 0.0821;  % resistance
L = 0.0308e-3; %inductance

V = 24;


% total inertia ([angular] acceleration terms)
J = J_m + J_g + J_s + (1/(N^2*n_g))*( 2*J_p + J_l + J_b + J_s ); 
A = J;

% Damping terms ([angular] velocity terms)
B = b_m + b + (b/(N^2*n_g));

% Force terms
% D = (r_p * mu_s * (m_l + m_b) * g * cos(angl) + r_p*(m_l + m_b) * g * sin(angl)) / N; 
D = (mu_s*cos(angl) + sin(angl)) * r_p*(m_l+m_b)*g/(N*n_g);

% UNCOMMENT NO LOAD CONDITION **
% A = J_m;
% B = b_m;
% D = 0;

 