% J.Williams
% University of Wisconsin-Madison
% EPD629: Powertrain Systems and Controls
% Assignement 10
% 2022-04-04

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assignment:
% This homework assignment is to create bode plots for an electronic
% throttle system and determine gains for a PI and a PID controller.
%  
% Part 1: Model Initialization File --- PI Control with Zero at -5rad/s (5 points)
%
% The remainder of the assignment is done via simulink.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note that this was assignment had a large copy-paste from the homework
% assignment instructions...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Close figures, clear variables and command window
clear;
close all;
clc;

% Assign symbolic variables
s = tf('s')  % Define s for symbolic transfer function format
% syms s L R J b K Kt Kb
  
% Calibration parameters
Kt = 0.063;  % Nm/A
Kb = Kt;     % V/(rad/s)
% L = 1e-3 was replaced by L = 1e-4 to keep poles real for this assignment
L = 1e-4;    % H  
R = 0.35;    % Ohm
J = 50e-6;   % kg*m^2
b = .001;    % Nm/(rad/s)
K = 0.1;     % Nm/rad
Tau = 0.001; % filter time constant
 
% Plant Model Transfer Function
G1 = Kt/(L*s+R)
G2 = 1/(J*s)
G2 = feedback(G2,b) % was G2 = G2/(1+G2*b) 
G2 = G2/s
G2 = feedback(G2,K) % was G2 = G2/(1+G2*K)
Gp = feedback(G1*G2,Kb*s) % was Gp = G1*G2/(1 + G1*G2*Kb*s)
                          % simplify(Gp)

% Plant poles
p = pole(Gp)
p1 = p(1) % Gp pole 1
p2 = p(2) % Gp pole 2
p3 = p(3) % Gp pole 3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Part 1: Proportional Control')
figure(1)
Gc_r1 = 1 % Gc for root locus with overall gain equal to one
rlocus(Gc_r1*Gp)
% This will be used to find the input for the gain used below...
Gain = 3.84  % Gain from breakaway point in root locus
Gc = Gain*Gc_r1 % Controller transfer function
sisotool(Gp, Gc) % Opens the Control System Designer tool

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Part 2: Integral Control')
figure(2)
Gc_r2 = 1/s % Gc for root locus with overall gain equal to one
rlocus(Gc_r2*Gp)
% This will be used to find the input for the gain used below...
Gain = 1.14  % Gain from breakaway point in root locus
Gc = Gain*Gc_r2 % Controller transfer function
sisotool(Gp, Gc) % Opens the Control System Designer tool

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Part 3: PI Control with Zero at -8.3691 rad/s')
% Note that -8.3691 should be p3
figure(3)
Gc_r3 = (s/(-p3)+1)/s % Gc for root locus with overall gain equal to one
rlocus(Gc_r3*Gp)
% This will be used to find the input for the gain used below...
Gain = 34.4  % Gain from breakaway point in root locus
Gc = Gain*Gc_r3 % Controller transfer function
sisotool(Gp, Gc) % Opens the Control System Designer tool

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Part 4: PID Control with Zeros at -256.9902 & -8.3691 rad/s')
% Note that -8.3691 should be p3, -256.9902 should be p2
figure(4)
Gc_r4 = (s/(-p2)+1)*(s/(-p3)+1)/(s*(Tau*s+1)) % Gc for root locus with overall gain equal to one
rlocus(Gc_r4*Gp)
% This will be used to find the input for the gain used below...
Gain = 119  % Gain from breakaway point in root locus
Gc = Gain*Gc_r4 % Controller transfer function
sisotool(Gp, Gc) % Opens the Control System Designer tool

%%%%% Need to do part 5 still


% Continue in Simulink...
