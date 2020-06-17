%%CoM
clear all
close all
clc

%Data
g = 9.81;
h_0 = 0.74;
w_0 = sqrt(9.81)/h_0;
A = [w_0,0;0,-w_0];
B = [-w_0;w_0];
C = eye(2);
D = [0;0];

LIPM = ss(A,B,C,D);
ns = 1;
deltaT = 0.01;
Tfinal = 1;
t = 0:deltaT:Tfinal;

%Time interval for steps
clock = 0.3;

alpha = 0.5;
x_c_0 = alpha*exp(-w_0*clock); %initial position of CoM
x_c_dot_0 = 0; %initial velocity of CoM
x_s_init = x_c_0 - x_c_dot_0/w_0;
x_u_init = x_c_0 + x_c_dot_0/w_0;

%Amplitudes of steps
% alpha = zeros(ns,1);
% 
% temp = exp(w_0*clock(1));
% alpha(1) = temp*x_u_init;
% for i=2:ns
%     alpha(i) = 0.5;
%     alpha(1) = alpha(1) - temp*alpha(i)*exp(-w_0*clock(i));
% end

%Control input
u = alpha*heaviside(t-clock);
%u = stairedSteps(deltaT,Tfinal,alpha,clock);

x_0 = [x_u_init;x_s_init];
y=lsim(LIPM,u,t,x_0);

%CoM trajectory
x_c = 1/2*(y(:,1)+y(:,2));
x_c_dot = w_0/2*(y(:,1)-y(:,2));

figure
plot(t,x_c,t,u,'r--','LineWidth',2)
xlabel('Time [s]')
ylabel('CoM Trajectory [m]')

figure
plot(t,x_c_dot,t,u,'r--','LineWidth',2)
xlabel('Time [s]')
ylabel('CoM Velocity [m/s]')

save('CoM_trajectory.mat','x_c','x_c_dot','deltaT','Tfinal');