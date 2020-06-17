%% Animation

load('animationInputs.mat');

%q1 = CoM x
%q2 = CoM y
%q3 = qtorso
%q4 = qfemStance
%q5 = qtibStance
%q6 = qfemSwing
%q7 = qtibSwing

anim = Animator.FiveLinkAnimator(t_log, q_log);
anim.pov = Animator.AnimatorPointOfView.West;
anim.Animate(true);
anim.isLooping = false;
anim.updateWorldPosition = true;
anim.endTime = 10;
conGUI = Animator.AnimatorControls();
conGUI.anim = anim;