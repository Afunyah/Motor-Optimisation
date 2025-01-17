% gear ratio
N = 1;  

% acceleration due to gravity
g = 9.81;

% masses
m_p = 0.025;    % pulley  1  
m_b = 0.05;    % belt 1.5
m_s = 0.025;    % shaft 0.6
m_l = 10;    % load 2

% radii
r_p = 0.5e-2;    % pulley
r_s = 0.3e-2;    % shaft


% inertia
J_m = 0.00000132;    % motor 
J_g = 0.0000;    % gearbox

J_p = m_p * r_p^2;  % pulley
J_l = m_l * r_p^2;  % load
J_b = m_b * r_p^2;  % belt
J_s = 0.5 * m_s * r_s^2;    % shaft

% Coefficients & Constants
b = 0e-6;      % shaft damping 1 * 10e-5
mu_v = 0.6e-6;   % viscous friction/damping
mu_s = 0.23;   % static friction (bed and belt) or couloumb?


Ke = 0.0168; % motor constant
R = 3.17;  % resistance
L = 0.333e-3; %inductance

V = 12;


% total inertia ([angular] acceleration terms)
J = J_m + J_g + J_s + (1/N^2)*( 2*J_p + J_l + J_b + J_s ); 
A = J;

% Damping terms ([angular] velocity terms)
B = b + (1/N^2)*(b + mu_v);

% Force terms
D = (r_p * mu_s * (m_l + m_b) * g) / N;  % I do not understand the role of static fricion ?or coulomb?
% D = 0; 


% UNCOMMENT NO LOAD CONDITION
% A = J_m;
% B = mu_v;
% D = 0;

