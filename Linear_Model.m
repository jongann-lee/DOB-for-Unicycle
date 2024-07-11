%% Parameter Setup
m_w = 2;
R = 0.15;
m = 3;
h = 0.3;
m_b = 1;
r_b = 0.15;
m_s = 1;
r_s = 0.15;
J_x = 0.1;
J_y = 0.1;
J_z = 0.02;
v_d = 3;

%% System state variables
%sigma1, sigma2, sigma3, sigma4, sigma5, sigma6, theta, gamma, psi, x_G, y_G, phi, alpha, beta;

%% Model Setup
% Notations
m_hat = m + m_b + m_s;
Q1 = 5*m_w*R^2 + 4*m_hat*(R+h)^2 + m_s*r_s^2 + 4*J_x;
Q2 = (3*m_w + 2*m_hat)*(m_b*r_b^2 + m_s*r_s^2 + 4 * J_y) + 12*m_hat*m_w*h^2;
Q3 = m_w*R^2 + m_b*r_b^2 + 4*J_z;

% A variables
A13 = (2*(3*m_w*R + 2*m_hat*(R + h))*v_d) / Q1;
A17 = (4*(m_w*R + m_hat*(R + h))*9.8) / Q1;
A28 = (-8*m_hat^2*h^2*9.8) / (R*Q2);
A31 = (-2*m_w*R*v_d) / Q3;
A48 = (4*(3*m_w + 2*m_hat)*m_hat*h*9.8) / Q2;

% B variables
B12 = (-4) / Q1;
B21 = (2*(4*m_hat*(R + h)*h + m_b*r_b^2 + m_s*r_s^2 + 4*J_y)) / (R^2*Q2);
B33 = (-4) / Q3;
B41 = (-4*(3*m_w*R + 2*m_hat*(R + h))) / (R*Q2);

% Linear Model
A_n = [0 0 A13 0 0 0 A17 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 A28 0 0 0 0 0 0;
    A31 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 A48 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    1 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 1 0 0 0 0 0 0 0 0 0 0;
    0 0 1 0 0 0 0 0 0 0 0 0 0 0;
    0 R 0 0 0 0 0 0 0 0 0 0 0 0;
    -R 0 0 0 0 0 0 0 v_d 0 0 0 0 0;
    0 1 0 0 0 0 0 0 0 0 0 0 0 0;
    -1 0 0 0 1 0 0 0 0 0 0 0 0 0;
    0 0 -1 0 0 1 0 0 0 0 0 0 0 0];
B_n = [0 B12 0;
    B21 0 0;
    0 0 B33;
    B41 0 0;
    0 2/(m_b*r_b^2) 0;
    0 0 2/(m_s*r_s^2);
    0 0 0;
    0 0 0;
    0 0 0;
    0 0 0;
    0 0 0;
    0 0 0;
    0 0 0;
    0 0 0];

C_n = zeros(3,14);
C_n(1,7) = 1; C_n(2,8) = 1; C_n(3,9) = 1;

%% Controller Design

Co = ctrb(A_n, B_n);
rank(Co) % system is not controllable
p_w = 15; d_w = 5;
p_b = 30; d_b = 10;
p_s = 15; d_s = 5;
% For the current PD controller
K = [0 0 0 d_w 0 0 0 p_w 0 0 0 0 0 0;
    d_b 0 0 0 0 0 p_b 0 0 0 0 0 0 0;
    0 0 d_s 0 0 0 0 0 p_s 0 0 0 0 0];

eigenvalues = eig(A_n+B_n*K);
%disp(eigenvalues)
%% DOB design

tau = 0.001;

a10 = 0.01; a20 = 0.01; a30 = 0.01;
a11 = 20; a21 = 20; a31 = 20;

A_a1 = [0 1;
    (-a10)/(tau^2) (-a11)/(tau)];
A_a2 = [0 1;
    (-a20)/(tau^2) (-a21)/(tau)];
A_a3 = [0 1;
    (-a30)/(tau^2) (-a31)/(tau)];

A_a = blkdiag(A_a1, A_a2, A_a3);
B_a = blkdiag((a10/(tau^2))*[0;1], (a20/(tau^2))*[0;1], (a30/(tau^2))*[0;1]);

%% Normal form stability
A_normal= [0 1 0 0 0 0;
    A17 0 0 0 0 A13;
    0 0 0 1 0 0;
    0 0 A48 0 0 0;
    0 0 0 0 0 1;
    0 A31 0 0 0 0];
B_normal = [0 0 0;
    0 B12 0;
    0 0 0;
    B41 0 0;
    0 0 0;
    0 0 B33];
K_normal = [0 0 p_w d_w 0 0;
    p_b d_b 0 0 0 0;
    0 0 0 0 p_s d_s];
eigenvalues_norm = eig(A_normal+B_normal*K_normal);
disp(eigenvalues_norm)



