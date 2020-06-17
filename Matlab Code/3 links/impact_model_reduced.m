%A MATLAB script to generate the equations of motion for a 5-link,
%planar, kneed, biped walker using the method of Lagrange, with 7
%gen coordinates.

%%%%%%%%%%%%%%
%IMPACT MODEL%
%%%%%%%%%%%%%%


clear all
close all

% -----------------------------------------------------------------
%
%  Model variables
%
% -----------------------------------------------------------------

% absolute joint angles and velocities
syms q_torso dq_torso real
syms q_fem1 dq_fem1 real
syms q_fem2 dq_fem2 real

% cartesian coordinates of the hip
syms pHx real 
syms pHy real

% linear velocities of the hip
syms vHx
syms vHy

% gravity
syms g real
g = double(subs(g,9.81));

% link lengths
syms L_torso real
syms L_fem real

L_torso = double(subs(L_torso,0.625)); %meters
L_fem = double(subs(L_fem,0.4));

% link masses
syms M_torso real 
syms M_fem real 

M_torso = double(subs(M_torso,20)); %kilos
M_fem = double(subs(M_fem,6.8));

% center of mass offsets (times mass??)
syms MY_torso real
syms MZ_torso real
syms MZ_fem real

MY_torso = double(subs(MY_torso,0)); %kilos times meter
MZ_torso = double(subs(MZ_torso,0.2));
MZ_fem = double(subs(MZ_fem,0.163));

% link inertias
syms XX_torso real
syms XX_fem real

XX_torso = double(subs(XX_torso,2.22)); %kilos times meters squared
XX_fem = double(subs(XX_fem,1.08));

% relative joint angles
syms q31L q32L q1L real
q  = [q31L; q32L; q1L];

% relative joint velocities
syms dq31L dq32L dq1L real
dq = [dq31L; dq32L; dq1L];


% -----------------------------------------------------------------
%
%  Change coordinates to absolute coordinates
%
% -----------------------------------------------------------------

T = [ 1   0   1;
      0   1   1;
      0   0   1];

% old absolute coordinates in terms of relative coords
q_new  = T*q;
dq_new = T*dq;

q_fem1  = q_new(1);
q_fem2  = q_new(2);
q_torso = q_new(3);
dq_fem1  = dq_new(1);
dq_fem2  = dq_new(2);
dq_torso = dq_new(3);

q_e = [q;pHx;pHy];
dq_e = [dq;vHx;vHy];
q_e_minus = q_e;
dq_e_minus = dq_e;

% -----------------------------------------------------------------
%
%  Calculate kinetic energy
%
% -----------------------------------------------------------------

% knee positions
p_hip = [pHx;pHy];
% p_knee1 =  p_hip + L_fem*[-sin(q_fem1); cos(q_fem1)];
% p_knee2 =  p_hip + L_fem*[-sin(q_fem2); cos(q_fem2)];

% hip and knee velocities
v_hip  = [vHx;vHy];
% v_knee1 = jacobian(p_knee1,q_e)*dq_e;
% v_knee2 = jacobian(p_knee2,q_e)*dq_e;

% relative angular velocties -- needed since the centers of masses
% of the links are one collocated with the link reference frames
R_torso = [cos(q_torso) -sin(q_torso);
	   sin(q_torso) cos(q_torso)];
v_torso = R_torso.'*v_hip*dq_torso;

R_fem1 = [cos(q_fem1) -sin(q_fem1);
	  sin(q_fem1) cos(q_fem1)];
v_fem1 = R_fem1.'*v_hip*dq_fem1;

R_fem2 = [cos(q_fem2) -sin(q_fem2);
	  sin(q_fem2) cos(q_fem2)];
v_fem2 = R_fem2.'*v_hip*dq_fem2;

% kinetic energy of links
KE_torso = 1/2*M_torso*(v_hip.'*v_hip) ...
	  + v_torso.'*[-MZ_torso; MY_torso]  ...
	  + 1/2*XX_torso*dq_torso^2;
KE_torso = simplify(KE_torso);
KE_fem1 = 1/2*M_fem*(v_hip.'*v_hip) ...
	 + v_fem1.'*[-MZ_fem; 0] ...
	 + 1/2*XX_fem*(dq_fem1)^2;
KE_fem1    = simplify(KE_fem1);
KE_fem2 = 1/2*M_fem*(v_hip.'*v_hip) ...
	 + v_fem2.'*[-MZ_fem; 0] ...
	 + 1/2*XX_fem*(dq_fem2)^2;
KE_fem2    = simplify(KE_fem2);

% total kinetic energy
KE = KE_torso + KE_fem1 + KE_fem2;
KE = simplify(KE);

% -----------------------------------------------------------------
%
%  Calculate potential energy
%
% -----------------------------------------------------------------

% positions of various members 
p_torso = p_hip ...
	  + 1/M_torso*[-sin(q_torso)*MZ_torso - cos(q_torso)*MY_torso;
		    cos(q_torso)*MZ_torso + sin(q_torso)*MY_torso];

p_fem1  = p_hip + MZ_fem/M_fem*[-sin(q_fem1); cos(q_fem1)];
p_fem2  = p_hip + MZ_fem/M_fem*[-sin(q_fem2); cos(q_fem2)];
p_foot1 =  p_hip + L_fem*[-sin(q_fem1); cos(q_fem1)];
p_foot2 =  p_hip + L_fem*[-sin(q_fem2); cos(q_fem2)];

% total potential energy
PE = g*(M_torso*p_torso(2) + M_fem*p_fem1(2) + M_fem*p_fem2(2));
PE = simplify(PE);

% -----------------------------------------------------------------
%
%  Calculate model matrices
%
% -----------------------------------------------------------------

% gravity vector
G_e = jacobian(PE,q_e).';
G_e = simplify(G_e);

% mass-inertial matrix
D_e = simplify(jacobian(KE,dq_e).');
D_e = simplify(jacobian(D_e,dq_e));

% Coriolis and centrifugal matrix
% syms C_e real
% n=max(size(q));
% for k=1:n
%   for j=1:n
%     C_e(k,j)=0*g;
%     for i=1:n
%       C_e(k,j)=C_e(k,j)+1/2*(diff(D_e(k,j),q(i)) + ...
% 			 diff(D_e(k,i),q(j)) - ...
% 			 diff(D_e(i,j),q(k)))*dq(i);
%     end
%   end
% end
% C_e=simplify(C_e);

% input matrix
% Phi_0 = [q_fem1-q_torso;
% 	 q_fem2-q_torso;
% 	 q_tib1-q_fem1;
% 	 q_tib2-q_fem2;
% 	 q_torso];
% B = jacobian(Phi_0,q);
% B = B.'*[eye(4,4);zeros(1,4)];

% swing foot force input matrix (F_ext = [F_T;F_N])
E = p_foot2;
dE = jacobian(E,q_e).';

%Relabeling matrix
R = [0 1 0;...
     1 0 0;...
     0 0 0];
 

Pi_g = inv([D_e, -dE;...
        dE', zeros(2)]);    
dq_e_plus = Pi_g(1:5,1:5)*D_e*dq_e_minus;


DELTA_dq = [R,zeros(3,2)]*Pi_g(1:5,1:5)*D_e;

%Switch relation
q_plus = R*q_e(1:3);
dq_plus = R*dq_e_plus(1:3);

