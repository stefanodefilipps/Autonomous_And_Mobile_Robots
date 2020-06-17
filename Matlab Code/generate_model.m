%A MATLAB script to generate the equations of motion for a 5-link,
%planar, kneed, biped walker using the method of Lagrange.
%
%Copyright (c) 2001 by Eric Westervelt and Jessy Grizzle.  This
%code may be freely used for noncommercial ends. If use of this
%code in part or in whole results in publication, proper citation
%must be included in that publication.  This code comes with no
%guarantees or support.
%

% This is for a robot with an upright torso, two legs and
% knees.  The model is for five degrees of freedom of the robot,
% with the stance foot (Foot1) pinned to the ground.  The notation
% used for the equation of motions are as in Robot Dynamics and
% Control by Spong and Vidyasagar (1989), page 142, Eq. (6.3.12)
%
%    D(q)ddq + C(q,dq)*dq + G(q) = B*tau + E2*F_external
%
% Note the following convention
%
%  fem = femur, tib = tibia, 1 = stance leg, 2 = swing leg
%
% For simplicity, the model is derived using absolute coorinates
% measured in the trigonometric sense from the vertical.  The model
% is trasnformed in to relative coordates plus one absolute
% coordinates via a coordinate transformation.
%
% Eric Westervelt, Mon Aug  6 17:08:49 EST 2001

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
syms q_tib1 dq_tib1 real
syms q_tib2 dq_tib2 real

% gravity
syms g real
g = double(subs(g,9.81));

% link lengths
syms L_torso real
syms L_fem real
syms L_tib real

L_torso = double(subs(L_torso,0.625)); %meters
L_fem = double(subs(L_fem,0.4));
L_tib = double(subs(L_tib,0.4));

% link masses
syms M_torso real 
syms M_fem real 
syms M_tib real

M_torso = double(subs(M_torso,20)); %kilos
M_fem = double(subs(M_fem,6.8));
M_tib = double(subs(M_tib,3.2));

% center of mass offsets (times mass??)
syms MY_torso real
syms MZ_torso real
syms MZ_fem real
syms MZ_tib real

MY_torso = double(subs(MY_torso,0)); %kilos times meter
MZ_torso = double(subs(MZ_torso,0.2));
MZ_fem = double(subs(MZ_fem,0.163));
MZ_tib = double(subs(MZ_tib,0.128));

% link inertias
syms XX_torso real
syms XX_fem real
syms XX_tib real

XX_torso = double(subs(XX_torso,2.22)); %kilos times meters squared
XX_fem = double(subs(XX_fem,1.08));
XX_tib = double(subs(XX_tib,0.93));

% relative joint angles
syms q31L q32L q41L q42L q1L real
q  = [q31L; q32L; q41L; q42L; q1L];

% relative joint velocities
syms dq31L dq32L dq41L dq42L dq1L real
dq = [dq31L; dq32L; dq41L; dq42L; dq1L];

% -----------------------------------------------------------------
%
%  Change coordinates to absolute coordinates
%
% -----------------------------------------------------------------

T = [ 1   0   0   0   1;
      0   1   0   0   1;
      1   0   1   0   1;
      0   1   0   1   1;
      0   0   0   0   1];

% old absolute coordinates in terms of relative coords
q_new  = T*q;
dq_new = T*dq;

q_fem1  = q_new(1);
q_fem2  = q_new(2);
q_tib1  = q_new(3);
q_tib2  = q_new(4);
q_torso = q_new(5);
dq_fem1  = dq_new(1);
dq_fem2  = dq_new(2);
dq_tib1  = dq_new(3);
dq_tib2  = dq_new(4);
dq_torso = dq_new(5);

% -----------------------------------------------------------------
%
%  Calculate kinetic energy
%
% -----------------------------------------------------------------

% hip and knee positions
p_hip  = L_fem*[sin(q_fem1); -cos(q_fem1)] ...
      + L_tib*[sin(q_tib1); -cos(q_tib1)];
p_knee1 = p_hip + L_fem*[-sin(q_fem1); cos(q_fem1)];
p_knee2 = p_hip + L_fem*[-sin(q_fem2); cos(q_fem2)];

% hip and knee velocities
v_hip  = jacobian(p_hip,q)*dq;
v_knee1 = jacobian(p_knee1,q)*dq;
v_knee2 = jacobian(p_knee2,q)*dq;

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

R_tib1 = [cos(q_tib1) -sin(q_tib1);
	  sin(q_tib1) cos(q_tib1)];
v_tib1 = R_tib1.'*v_knee1*dq_tib1;

R_tib2 = [cos(q_tib2) -sin(q_tib2);
	  sin(q_tib2) cos(q_tib2)];
v_tib2 = R_tib2.'*v_knee2*dq_tib2;

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
KE_tib1 = 1/2*M_tib*(v_knee1.'*v_knee1) ...
	 + v_tib1.'*[-MZ_tib; 0] ...
	 + 1/2*XX_tib*(dq_tib1)^2;
KE_tib1 = simplify(KE_tib1);
KE_tib2 = 1/2*M_tib*(v_knee2.'*v_knee2) ...
	 + v_tib2.'*[-MZ_tib;0] ...
	 + 1/2*XX_tib*(dq_tib2)^2;
KE_tib2 = simplify(KE_tib2);

% total kinetic energy
KE = KE_torso + KE_fem1 + KE_fem2 + KE_tib1 + KE_tib2;
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
p_tib1  = p_knee1 + MZ_tib/M_tib*[-sin(q_tib1); cos(q_tib1)];
p_tib2  = p_knee2 + MZ_tib/M_tib*[-sin(q_tib2); cos(q_tib2)];
p_foot1 = p_knee1 + L_tib*[-sin(q_tib1); cos(q_tib1)];
p_foot2 = p_knee2 + L_tib*[-sin(q_tib2); cos(q_tib2)];

% total potential energy
PE = g*(M_torso*p_torso(2) + M_fem*p_fem1(2) + M_fem*p_fem2(2) + ...
	M_tib*p_tib1(2) + M_tib*p_tib2(2));
PE = simplify(PE);

% -----------------------------------------------------------------
%
%  Calculate model matrices
%
% -----------------------------------------------------------------

% gravity vector
G = jacobian(PE,q).';
G = simplify(G);

% mass-inertial matrix
D = simplify(jacobian(KE,dq).');
D = simplify(jacobian(D,dq));

% Coriolis and centrifugal matrix
syms C real
n=max(size(q));
for k=1:n
  for j=1:n
    C(k,j)=0*g;
    for i=1:n
      C(k,j)=C(k,j)+1/2*(diff(D(k,j),q(i)) + ...
			 diff(D(k,i),q(j)) - ...
			 diff(D(i,j),q(k)))*dq(i);
    end
  end
end
C=simplify(C);

% input matrix
Phi_0 = [q_fem1-q_torso;
	 q_fem2-q_torso;
	 q_tib1-q_fem1;
	 q_tib2-q_fem2;
	 q_torso];
B = jacobian(Phi_0,q);
B = B.'*[eye(4,4);zeros(1,4)];

% swing foot force input matrix (F_ext = [F_T;F_N])
E = [p_foot2];
dE = jacobian(E,q).';

%Drift f
f1 = dq;
f2 = D\(-C*dq-G);
f = [f1;f2];

%Input vector field
g1 = zeros(5,4);
g2 = D\B;
g3 = [g1;g2];
