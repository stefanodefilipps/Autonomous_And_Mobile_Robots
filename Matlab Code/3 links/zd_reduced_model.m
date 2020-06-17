%%Compute Zero Dynamics%%

clear all
close all

load('model_reduced.mat');

%Coordinate transformation
syms z1 z2 z3 z4 real
syms x1 real
syms x2 real
z = [z1;z2;z3;z4];
x = [x1;x2];
new_coord = [z;x];

%9th coord transf
%c = [-1,0,-1/2,0,-1];
%given the current structure of c when whe have 5 generalized coordinates i
%would do

c = [-1,-1/2,-1];

theta = c*q;

%10th coord transf
gamma = D(3,:)*dq;

%Lie derivatives in the original coord
L_f_theta = [c,zeros(1,3)]*f;
L_f_gamma = -G(3,:);

%q_0 - Initial pose
eqn1 = p_foot2(2) == 0;
eqn2 = p_foot1(2) == 0;
eqn3 = p_foot2(1) == 0.2;
[q_0_1, q_0_2, q_0_3] = vpasolve([eqn1,eqn2,eqn3],q,...
    [160*pi/180,200*pi/180,-20*pi/180]);

q_0 = [single(q_0_1);...
    single(q_0_2);...
    single(q_0_3)];

q_0_degree = q_0*180/3.14;

%OUTPUT

%Define h_0
H_0 = [eye(2),zeros(2,1)];
h_0 = H_0*q;
H = [H_0', c'];
theta_minus = c*q_0;
M = 6;
R = [0 1 0;...
     1 0 0;...
     0 0 0];
theta_plus = c*R*q_0;

%Define h_d composed s(q)
alpha_matrix = sym('alpha_matrix',[4,7],'real');
h_d = zeros(2,1); h_d = sym(h_d);
s=(x1-theta_plus)/(theta_minus-theta_plus);
for i=1:2
    for j=1:(M+1)
        h_d(i) = h_d(i) + alpha_matrix(i,j)*...
            factorial(M)/(factorial(j)*factorial(M-j+1))*...
            s^(j-1)*(1-s)^(M-j+1);     
    end
end
h_d_q = subs(h_d,x1,c*q);

%Final output
h = h_0 - h_d_q;

dh = jacobian(h,q);


% alpha_0_transform = H*R/H*[alpha_matrix(:,5);theta_minus];
% alpha_0_transform = alpha_0_transform(1:4);
% gamma_0 = D(5,:);
% gamma_0 = simplify(gamma_0); 
% foo_matrix = [dh;gamma_0];
% foo_matrix_0 = subs(foo_matrix,q,q_0);
% q_dot_0_minus = foo_matrix_0\[zeros(4,1);1]; 
% q_dot_plus = subs(DELTA_dq,q,q_0)*[q_dot_0_minus;0;0];
% q_dot_plus = q_dot_plus(1:5);
% alpha_1_transform = (theta_minus - theta_plus)/(M*c*...
%     q_dot_plus)*H_0*q_dot_plus +...
%     alpha_0_transform;
% h = subs(h,alpha_matrix(:,2),[alpha_1_transform]);
% h = subs(h,alpha_matrix(:,1),[alpha_0_transform]);

% %Inverse diffeomorphism
x2q = H\[h_d;x1];
gamma_0 = D(3,:);
foo_matrix = [dh;gamma_0];
q2dq = (foo_matrix)\[zeros(2,1);1];
q2dq = q2dq*x2;
x2dq = subs(q2dq,q,x2q);
% kappa1 = c*q2dq;
% 
% % % %Final zero dynamics in x coordinates
L_f_theta_X = subs(L_f_theta,dq,x2dq);
L_f_gamma_X = subs(L_f_gamma,q,x2q);


%% Feedback linearization

state = [q;dq];

L_f_h = jacobian(h,state)*f;
L_f_squared_h = jacobian(L_f_h,state)*f;


%Original coordinates q
L_f_squared_h_X = subs(L_f_squared_h,state,[x2q;x2dq]);
% 
% %ZD coordinates
% L_f_squared_h_X = vpa(L_f_squared_h_X,3);
% 
L_g_L_f_h = jacobian(L_f_h,state)*g3;
% 
% %Original coordinates q
% L_g_L_f_h = vpa(L_g_L_f_h,3);
% L_g_L_f_h_X = subs(L_g_L_f_h,q,x2q);
% 
% %ZD coordintes
% L_g_L_f_h_X = vpa(L_g_L_f_h_X,3);

%Control Input
% u = -L_g_L_f_h\L_f_squared_h;
% I noticed that I need to compute it this way otherwise if I use \ it is
% stuck at computing
L_g_L_f_h_inv = inv(L_g_L_f_h);
u = -L_g_L_f_h_inv*L_f_squared_h;