%% Define parameter Symbols
syms m_w R m h J_x J_y J_z m_b r_b m_s r_s g M_w M_b M_s delta...
    x_G y_G z_G psi theta phi gamma alpha beta x_dot_G y_dot_G z_dot_G...
    psi_dot theta_dot phi_dot gamma_dot alpha_dot beta_dot...
    sigma1 sigma2 sigma3 sigma4 sigma5 sigma6...
    sigma_dot_1 sigma_dot_2 sigma_dot_3 sigma_dot_4 sigma_dot_5 sigma_dot_6;

%% Frame Transformation
T_02 = [cos(psi), sin(psi), 0;...
    -cos(theta)*sin(phi), cos(psi)*cos(theta), sin(theta);...
    sin(psi)*sin(theta), -cos(psi)*sin(theta),cos(theta)];
T_23=[cos(gamma), 0, -sin(gamma);...
    0,1,0;...
    sin(gamma),0,cos(gamma)];
T_03=T_23*T_02;

%% Constraints derivation
vG_F0=[x_dot_G; y_dot_G; z_dot_G];
omega_w_F2=T_02*[0;0;psi_dot] + [theta_dot; phi_dot;0];
rGP_F2=[0;0;-R];
vP_F0=vG_F0+inv(T_02)*cross(omega_w_F2,rGP_F2); 
VP_F0=simplify(vP_F0); 
constraints=vP_F0==[0;0;0];

% kinematic constraints
KC=[vP_F0(1); vP_F0(2)]==[0;0];
GC=z_G==R*cos(theta);

%% Pseudo Velocities
omega_F3=T_03*[0;0;psi_dot]+T_23*[theta_dot; 0; 0]+[0; gamma_dot;0];
omega_b_F3=omega_F3+[alpha_dot;0;0];
omega_s_F3=omega_F3 + [0;0; beta_dot];

% pseudo velocities
pseudo_vel=[sigma1; sigma2; sigma3; sigma4; sigma5; sigma6]==...
[omega_w_F2; omega_F3(2); omega_b_F3(1); omega_s_F3(3)];

%% Kinematics Equations
A=[1,0,-R*cos(psi) *sin(theta), -R*cos(theta) *sin(psi), -R*cos(psi),0,0,0;... 
    0,1,-R*sin(psi) *sin(theta), R*cos(psi) *cos(theta), -R*sin(psi),0,0,0;... 
    0,0,0,1,0,0,0,0;...
    0,0,sin(theta),0,1,0,0,0;...
    0,0,cos(theta),0,0,0,0,0;... 
    0,0,sin(theta),0,0,1,0,0;...
    0,0,-cos(theta)*sin(gamma), cos(gamma),0,0,1,0;...
    0,0, cos(gamma)*cos(theta), sin(gamma),0,0,0,1];
det(A);
B=[x_dot_G;y_dot_G;psi_dot; theta_dot; phi_dot; gamma_dot; alpha_dot; beta_dot]; 
C=[0;0; sigma1; sigma2; sigma3; sigma4; sigma5; sigma6]; 
A*B==C;
B_soln=inv(A)*C;
% kinematic equations 
EOM_p1=B==B_soln;

%% Acceleration Energy of Wheel 
vG_F2=[0;0;0]-cross (omega_w_F2,rGP_F2);
vG_F2=subs(vG_F2,B,B_soln);
omega_20_F2=T_02* [0;0;psi_dot]+[theta_dot;0;0];
omega_20_F2=subs(omega_20_F2,B,B_soln);
aG_F2=[R*sigma_dot_2; -R*sigma_dot_1;0]+cross (omega_20_F2, vG_F2);
S_wheel_T1=(1/2)*m_w*(aG_F2 (1)^2+aG_F2 (2)^2+aG_F2(3)^2);

omega_w_F2=[sigma1; sigma2; sigma3];
alpha_w_F2=[sigma_dot_1; sigma_dot_2; sigma_dot_3]+cross(omega_20_F2, omega_w_F2);
J_w_F2=[(1/4)*m_w*R*2,0,0;...
    0, (1/2)*m_w*R*2,0;...
    0,0, (1/4)*m_w*R^2];
S_wheel_T2=(1/2)* [alpha_w_F2(1), alpha_w_F2(2), alpha_w_F2(3)]*J_w_F2*alpha_w_F2;

H_w_F2=J_w_F2*omega_w_F2;
S_wheel_T3=det([alpha_w_F2(1), alpha_w_F2(2), alpha_w_F2(3); ...
    omega_w_F2(1), omega_w_F2(2), omega_w_F2(3); ...
    H_w_F2(1),H_w_F2(2),H_w_F2(3)]);
S_wheel=S_wheel_T1+S_wheel_T2+S_wheel_T3;

%% Acceleration Energy of Body
vG_F3=T_23*vG_F2;
omega_F3=subs(omega_F3,B,B_soln);
rGB_F3=[0;0;h];
vB_F3=vG_F3+cross (omega_F3,rGB_F3);
vB_F3_local_der(1,1)=h*sigma_dot_4+R*(sigma_dot_2*cos (gamma)-sigma2*gamma_dot*sin(gamma));
vB_F3_local_der(2,1)=-R*sigma_dot_1-h*(sigma_dot_1*cos(gamma)...
    -sigma1*gamma_dot*sin(gamma)-sigma_dot_3*sin(gamma)... 
    -sigma3*gamma_dot*cos (gamma));
vB_F3_local_der(3,1) = R*(sigma_dot_2*sin(gamma)+sigma2*gamma_dot*cos(gamma)); 
vB_F3_local_der=subs(vB_F3_local_der,B,B_soln);
omega_30_F3=omega_F3;
aB_F3=vB_F3_local_der+cross(omega_30_F3, vB_F3);
S_body_T1= (1/2)*m* (aB_F3 (1)^2+aB_F3 (2)^2+aB_F3(3)^2);

omega_F3_local_der(1,1)=sigma_dot_1*cos(gamma)-sigma1*gamma_dot*sin(gamma)...
    -sigma_dot_3*sin(gamma)-sigma3*gamma_dot*cos(gamma);

omega_F3_local_der(2,1)=sigma_dot_4; 
omega_F3_local_der(3,1)=sigma_dot_3*cos(gamma)-sigma3*gamma_dot*sin(gamma)...
    +sigma_dot_1*sin(gamma)+sigma1*gamma_dot*cos(gamma); 
omega_F3_local_der=subs(omega_F3_local_der,B,B_soln); 
alpha_F3=omega_F3_local_der+cross(omega_30_F3, omega_F3);
J_F3=[J_x,0,0; 0,J_y,0;0,0,J_z];
S_body_T2=(1/2)*[alpha_F3(1), alpha_F3(2), alpha_F3(3)]*J_F3*alpha_F3;

H_F3=J_F3*omega_F3;
S_body_T3=det([alpha_F3(1), alpha_F3(2), alpha_F3(3); ...
    omega_F3(1), omega_F3(2), omega_F3(3);... 
    H_F3(1),H_F3(2),H_F3(3)]);
S_body=S_body_T1+S_body_T2+S_body_T3;

%% Acceleration Energy of Balancing Flywheel 
S_balance_T1=(1/2)*m_b*(aB_F3 (1)^2+aB_F3(2)^2+aB_F3 (3)^2);

omega_b_F3=subs(omega_b_F3,B,B_soln);
omega_b_F3_local_der(1,1)=sigma_dot_5;
omega_b_F3_local_der(2,1)=sigma_dot_4;
omega_b_F3_local_der(3,1)=sigma_dot_3*cos (gamma)... 
    -sigma3*gamma_dot*sin(gamma) +sigma_dot_1*sin(gamma)... 
    +sigma1*gamma_dot*cos(gamma);
omega_b_F3_local_der=subs(omega_b_F3_local_der,B,B_soln); 
alpha_b_F3=omega_b_F3_local_der+cross (omega_30_F3, omega_b_F3);
J_b_F3 = [(1/2)*m_b*r_b^2,0,0;0, (1/4)*m_b*r_b^2,0;0,0, (1/4)*m_b*r_b^2]; 
S_balance_T2=(1/2)*[alpha_b_F3(1), alpha_b_F3(2), alpha_b_F3(3)] * J_b_F3*alpha_b_F3;

H_b_F3=J_b_F3*omega_b_F3;
S_balance_T3=det([alpha_b_F3(1), alpha_b_F3(2), alpha_b_F3(3);...
    omega_b_F3(1), omega_b_F3(2), omega_b_F3(3); ...
    H_b_F3(1),H_b_F3(2),H_b_F3(3)]);

S_balance=S_balance_T1+S_balance_T2+S_balance_T3;

%% Acceleration Energy of Steering Flywheel 
S_steer_T1=(1/2)*m_s* (aB_F3 (1)^2+aB_F3(2)^2+aB_F3(3)^2);
omega_s_F3=subs(omega_s_F3,B,B_soln);
omega_s_F3_local_der(1,1)=sigma_dot_1*cos (gamma)... 
    -sigma1*gamma_dot*sin(gamma)-sigma_dot_3*sin(gamma)... 
    -sigma3*gamma_dot*cos(gamma); 
omega_s_F3_local_der(2, 1)=sigma_dot_4; 
omega_s_F3_local_der(3, 1)=sigma_dot_6;
omega_s_F3_local_der=subs(omega_s_F3_local_der,B,B_soln); 
alpha_s_F3=omega_s_F3_local_der+cross(omega_30_F3, omega_s_F3); 
J_s_F3= [(1/4)*m_s*r_s^2,0,0;0, (1/4) *m_s*r_s^2,0;0,0, (1/2)*m_s*r_s^2]; 
S_steer_T2=(1/2)*[alpha_s_F3(1), alpha_s_F3(2), alpha_s_F3(3)]...
    *J_s_F3*alpha_s_F3;

H_s_F3=J_s_F3*omega_s_F3;
S_steer_T3=det([alpha_s_F3(1), alpha_s_F3(2), alpha_s_F3(3); ... 
    omega_s_F3(1), omega_s_F3(2), omega_s_F3(3); ... 
    H_s_F3(1),H_s_F3(2),H_s_F3(3)]);
S_steer=S_steer_T1+S_steer_T2+S_steer_T3;

%% Total Acceleration Energy 
S=S_wheel+S_body+S_balance+S_steer;

%% Virtual Power
% virtual power of rolling wheel
G_w_F2=T_02*[0;0;-m_w*g];
M_w_F2=[0; M_w; 0];
M_w_F3=T_23*M_w_F2;
delta_P_wheel=[G_w_F2(1),G_w_F2(2), G_w_F2(3)]...
    *[delta*vG_F2(1); delta*vG_F2(2); delta*vG_F2(3)]+... 
    [M_w_F2(1),M_w_F2(2),M_w_F2(3)]...
    *[delta*omega_w_F2(1); delta*omega_w_F2(2); delta*omega_w_F2(3)];

% virtual power of balancing flywheel
G_b_F3=T_03*[0;0;-m_b*g];
M_b_F3= [M_b; 0; 0];
delta_P_balance=[G_b_F3(1), G_b_F3(2), G_b_F3(3)]...
    * [delta*vB_F3(1); delta*vB_F3(2) ; delta*vB_F3(3)]+... 
    [M_b_F3(1),M_b_F3(2), M_b_F3(3)]...
    * [delta*omega_b_F3(1); delta*omega_b_F3(2) ; delta*omega_b_F3(3)];

% virtual power of steering flywheel 
G_s_F3=T_03* [0;0;-m_s*g];
M_s_F3=[0;0;M_s];
delta_P_steer=[G_s_F3(1), G_s_F3(2), G_s_F3(3)]...
    *[delta*vB_F3(1); delta*vB_F3(2); delta*vB_F3(3)]+...
    [M_s_F3(1), M_s_F3(2), M_s_F3(3)]...
    * [delta*omega_s_F3(1); delta*omega_s_F3(2); delta*omega_s_F3(3)];

% virtual power of body
G_F3=T_03* [0;0;-m*g];
delta_P_body=[G_F3(1), G_F3(2), G_F3(3)]...
    * [delta*vB_F3(1); delta*vB_F3(2) ; delta*vB_F3(3)]+...
    [-M_w_F3(1), -M_w_F3(2), -M_w_F3(3)]...
    * [delta*omega_F3(1); delta*omega_F3(2); delta* omega_F3(3)]+... 
    [-M_b_F3(1), -M_b_F3(2), -M_b_F3(3)]...
    *[delta*omega_F3(1); delta*omega_F3(2); delta*omega_F3(3)]+... 
    [-M_s_F3(1), -M_s_F3(2), -M_s_F3(3)]...
    *[delta*omega_F3(1); delta*omega_F3(2); delta*omega_F3(3)];

% total virtual power
delta_P=delta_P_wheel+delta_P_balance+delta_P_steer+delta_P_body;

%% Pseudo forces
Pi_1=R*g* (m+m_w+m_b+m_s)*sin(theta)...
    +g*h* (m+m_b+m_s) *cos (gamma) *sin(theta)-M_b*cos (gamma)-M_s*sin(gamma); 
Pi_2=M_w;
Pi_3=M_b*sin(gamma)-M_s*cos (gamma)-g*h*(m+m_b+m_s)*sin(gamma)*sin(theta);
Pi_4=g*h* (m+m_b+m_s)*cos(theta) *sin(gamma)-M_w;
Pi_5=M_b;
Pi_6=M_s;
%% Appell Equations
Appell_1=simplify(diff(S, sigma_dot_1)) == Pi_1; 
Appell_2=simplify(diff(S, sigma_dot_2)) == Pi_2;
Appell_3=simplify(diff(S, sigma_dot_3)) == Pi_3; 
Appell_4=simplify(diff(S, sigma_dot_4)) == Pi_4; 
Appell_5=simplify(diff(S, sigma_dot_5)) == Pi_5; 
Appell_6=simplify(diff(S, sigma_dot_6)) == Pi_6;
[sigma_dot_1_soln, sigma_dot_2_soln, sigma_dot_3_soln,...
    sigma_dot_4_soln, sigma_dot_5_soln, sigma_dot_6_soln]=...
    solve([Appell_1, Appell_2, Appell_3, Appell_4, Appell_5, Appell_6],... 
    [sigma_dot_1, sigma_dot_2, sigma_dot_3, sigma_dot_4, sigma_dot_5, sigma_dot_6]);
%% EOM
% Appell equations
f_1 = simplify(expand(sigma_dot_1_soln), 'Steps',50);
f_2 = simplify(expand(sigma_dot_2_soln), 'Steps',50);
f_3 = simplify(expand(sigma_dot_3_soln), 'Steps',50);
f_4 = simplify(expand(sigma_dot_4_soln), 'Steps',50);
f_5 = simplify(expand(sigma_dot_5_soln), 'Steps',50);
f_6 = simplify(expand(sigma_dot_6_soln), 'Steps',50);
EOM_p2=...
[sigma_dot_1==f_1; ... % f_1 
sigma_dot_2==f_2; ... % f_2 
sigma_dot_3==f_3; ... %f_3
sigma_dot_4==f_4; ... %f_4
sigma_dot_5==f_5; ... %f_5
sigma_dot_6==f_6];   %f_6
% EOM
EOM=[EOM_p2;EOM_p1];

%% Find G
[num1, den1] = numden(f_1);
G11 = 0;
num1_b = diff(num1, M_b);
G12 = num1_b/den1;
num1_s = diff(num1,M_s);
G13 = num1_s/den1;

temp2 = f_4 - f_3*tan(theta) - sigma3/((cos(theta))^2);
%dot2 = simplify(expand(temp2), 'Steps',50);
[num2,den2] = numden(temp2);
den2 = simplify(expand(den2), 'Steps',50);
num2_w = diff(num2,M_w);
num2_w = simplify(expand(num2_w), 'Steps',50);
G21 = num2_w/den2;
num2_b = diff(num2,M_b);
num2_b = simplify(expand(num2_b), 'Steps',50);
G22 = num2_b/den2;
num2_s = diff(num2,M_s);
num2_s = simplify(expand(num2_s), 'Steps',50);
G23 = num2_s/den2;

temp3 = f_3/(cos(theta)) + sigma3*tan(theta);
dot3 = simplify(expand(temp3), 'Steps',50);
[num3,den3] = numden(dot3);
G31 = 0;
num3_b = diff(num3, M_b);
num3_b = simplify(expand(num3_b), 'Steps',50);
G32 = num3_b/den3; 
num3_s = diff(num3,M_s);
num3_s = simplify(expand(num3_s), 'Steps',50);
G33 = num3_s/den3;

G = [G11 G12 G13;
    G21 G22 G23;
    G31 G32 G33];

%% Assumption 3 solver

syms state1 state2 state3;
state = [state1; state2; state3];
K = eye(3);
G_minus = diag([0.1 0.1 0.1]);
G_plus = diag([10 10 10]);
BigPi = inv((G_plus + G_minus)/2);

MI = (G*K*state - G_minus*state)'*BigPi^2*(G*K*state - G_plus*state) <= 0;

constraints = [theta > -1.3, theta < 1.3, ...
    gamma > -1.3, gamma < 1.3,...
    psi > -1.3, psi < 1.3];

sol3 = solve(MI, constraints);

