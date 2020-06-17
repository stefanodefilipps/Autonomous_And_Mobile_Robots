%% Kinematic tracking
clear all
close all
clc

%Log inputs for reproducing animation
load('animationInputs.mat');

%Inputs for CoM trajectory
load('CoM_trajectory.mat');

%Direct Kinematics for the hip
load('p_hip.mat');

%Direct Kinematics for the swing foot
load('p_swing.mat');


%% Reference for the hip angle
q_torso_ref = -q_log_R(3,:);
t_torso = 0:Tfinal/20:Tfinal;
t_torso_new = 0:Tfinal/(20.2*5):Tfinal;
q_torso_ref = interp1(t_torso,q_torso_ref,t_torso_new,'spline');
q_log_degree = 180/pi*q_log(3:7,:);
q_torso_dot_ref = zeros(1,101);
Tf = size(x_c,1);

for i=1:Tf
    q_torso_dot_ref(i) = (q_torso_ref(i+1)-q_torso_ref(i))/(Tfinal/(20.2*5));
end

%% Reference for the swing foot

%x coordinate
%The swing foot should start from -0.37 m
%and land at 0.4 m
x_sw_0 = -0.37;
x_sw_F = 0.37;
x_sw_ref = x_sw_0 + t_torso_new*(x_sw_F-x_sw_0);
x_sw_dot_ref = zeros(1,101);
for i=1:Tf
    x_sw_dot_ref(i) = (x_sw_ref(i+1)-x_sw_ref(i))/(Tfinal/(20.2*5));
end

%z coordinate
%The swing foot maximum height should be 0.1 m
z_sw_MAX = 0.1;
z_sw_ref = -4*z_sw_MAX/Tfinal^2*t_torso_new.*...
    (t_torso_new-Tfinal);
z_sw_dot_ref = zeros(1,101);
for i=1:Tf
    z_sw_dot_ref(i) = (z_sw_ref(i+1)-z_sw_ref(i))/(Tfinal/(20.2*5));
end

swing_ref = [x_sw_ref;z_sw_ref];
swing_vel = [x_sw_dot_ref;z_sw_dot_ref];



%% Tracking

%Initial pose
q0 = [-0.1718;2*pi-2.42;-0.6430;2*pi-2.9432;-0.5858];

p_swing = p_foot2;

q_sym = symvar(p_swing);
q_sym = q_sym([1:2,4,3,5]);
p_hip_0 = double(subs(p_hip,q_sym,q0'));
p_swing_0 = double(subs(p_swing,q_sym,q0'));
delta_com_traj = abs(x_c(1)-p_hip_0(1));
x_c = x_c - delta_com_traj;

dir_kin = [p_hip;q_sym(1);p_swing];

jac = jacobian(dir_kin,q_sym);

z_c = p_hip_0(2)*ones(Tf,1);
z_c_dot = zeros(Tf,1);
q_dot = zeros(5,Tf);
q = zeros(5,Tf+1);
q(:,1) = q0;
q_actual = q0;
gain = 0;
com_ref = [x_c';z_c'];
stack_pos = [com_ref;q_torso_ref(1:end-1);swing_ref(:,1:end-1)];
com_vel = [x_c_dot';z_c_dot'];
stack_vel = [com_vel;q_torso_dot_ref;swing_vel];

stack_num = zeros(5,Tf+1);

for i=1:Tf
    j_num = single(subs(jac,q_sym,q_actual'));
    stack_num(:,i) = single(subs(dir_kin,q_sym,q_actual'));
    q_dot(:,i) = j_num\(stack_vel(:,i)+gain*(stack_pos(:,i) - ...
        stack_num(:,i)));
    q(:,i+1) = q(:,i) + deltaT*q_dot(:,i);
    q_actual = q(:,i+1);
end

%% Plots

%Positions

t = 0:deltaT:Tfinal;
figure
for i=1:5
    subplot(3,2,i)
    plot(t,q(i,1:end-1));
    xlabel('time [s]');
    ylabelstring = sprintf('q_%i [rad]',i);
    ylabel(ylabelstring);
    title("Joint Angles");
    grid
end

figure
for i=1:5
    subplot(3,2,i)
    plot(t,stack_num(i,1:end-1),t,stack_pos(i,:),'r--');
    xlabel('time [s]');
    ylabelstring = sprintf('[rad OR m]',i);
    ylabel(ylabelstring);
    title("Positions");
    grid
end


%Velocities

stack_num_vel = zeros(5,Tf);
for i=1:Tf
    stack_num_vel(:,i) = (stack_num(:,i+1)-stack_num(:,i))/(Tfinal/(20.2*5));
end

figure
for i=1:5
    subplot(3,2,i)
    plot(t,q_dot(i,:));
    xlabel('time [s]');
    ylabelstring = sprintf('q_{dot}_%i [rad/s]',i);
    ylabel(ylabelstring);
    title("Joint Velocities");
    grid
end

figure
for i=1:5
    subplot(3,2,i)
    plot(t(1:end-1),stack_num_vel(i,1:end-1),...
        t(1:end-1),stack_vel(i,1:end-1),'r--');
    xlabel('time [s]');
    ylabelstring = sprintf('[rad/s OR m/s]',i);
    ylabel(ylabelstring);
    title("Velocities");
    grid
end

%% New simulation
%Changing coordinates for simulation
new_q = q;
new_q(1,:) = -q(1,:);
new_q(3,:) = -q(3,:); new_q(5,:) = -q(5,:);
new_q(2,:) = -q(2,:) + 2*pi;
new_q(4,:) = -q(4,:) + 2*pi;

%Number of steps
num_steps = 10;

q_LOG_0 = [stack_num(1:2,1:end-1);new_q(:,1:end-1)];

time_vec = zeros(1,Tf,num_steps); 
time_vec(:,:,1) = t;
q_final = zeros(7,Tf,num_steps);

for i=1:num_steps
    if i==1
        q_temp = q_LOG_0;
        q_final(:,:,i) = q_temp;
    else
       if mod(i,2)==0
            q_temp = q_final([1:3,6:7,4:5],:,i-1); % symmetric leg (switching)
            q_temp(1,:) = q_temp(1,:) + ...
                repmat((q_final(1,end,i-1)-...
                q_final(1,1,i-1)),1,Tf);
            q_final(:,:,i) = q_temp;
       else
            q_temp = q_final([1:3,6:7,4:5],:,i-1);
            q_temp(1,:) = q_temp(1,:) + ...
                repmat((q_final(1,end,i-1)-...
                q_final(1,1,i-1)),1,Tf);
            q_final(:,:,i) = q_temp;
       end
       time_vec(:,:,i) = time_vec(:,:,i-1) + time_vec(:,end,1);
    end
end

t_log_tot = time_vec(:,:,1);
q_LOG_tot = q_final(:,:,1);

if num_steps > 1
    for i=2:num_steps
        t_log_tot = [t_log_tot,time_vec(:,:,i)];
        q_LOG_tot = [q_LOG_tot,q_final(:,:,i)];
    end
end

anim = Animator.FiveLinkAnimator(t_log_tot, q_LOG_tot);
anim.pov = Animator.AnimatorPointOfView.West;
anim.Animate(true);
anim.isLooping = false;
anim.updateWorldPosition = true;
anim.endTime = Tfinal*num_steps;
conGUI = Animator.AnimatorControls();
conGUI.anim = anim;

