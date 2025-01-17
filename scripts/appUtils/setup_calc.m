g = 9.81;
angl = deg2rad(angl);

J_g = 0.0000;

J_p = m_p * r_p^2;  % pulley
J_l = m_l * r_p^2;  % load
J_b = m_b * r_p^2;  % belt
J_s = 0.5 * m_s * r_s^2;    % shaft

% total inertia ([angular] acceleration terms)
J = J_m + J_g + J_s + (1/(N^2*n_g))*( 2*J_p + J_l + J_b + J_s ); 
A = J;

% Damping terms ([angular] velocity terms)
B = b_m + b + (b/(N^2*n_g));

% Force terms
% D = (r_p * mu_s * (m_l + m_b) * g * cos(angl) + r_p*(m_l + m_b) * g * sin(angl)) / N; 
D = (mu_s*cos(angl) + sin(angl)) * r_p*(m_l+m_b)*g/(N*n_g);
